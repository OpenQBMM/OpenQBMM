/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "univariatePopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(univariatePopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        univariatePopulationBalance, 
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceModels::univariatePopulationBalance
::univariatePopulationBalance
(
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    populationBalanceModel(dict, U, phi),
    quadrature_(U.mesh()),
    aggregation_(dict.lookup("aggregation")),
    breakup_(dict.lookup("breakup")),
    aggregationKernel_
    (
        Foam::populationBalanceSubModels::aggregationKernel::New
        (
            dict.subDict("aggregationKernel")
        )
    ),
    breakupKernel_
    (
        Foam::populationBalanceSubModels::breakupKernel::New
        (
            dict.subDict("breakupKernel")
        )
    ),
    daughterDistribution_
    (
        Foam::populationBalanceSubModels::daughterDistribution::New
        (
            dict.subDict("daughterDistribution")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceModels::univariatePopulationBalance
::~univariatePopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::populationBalanceModels::univariatePopulationBalance::advectMoments()
{
    Info<< "\nAdvecting moments." << endl;

    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("nei", dimless, -1.0)
    );

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("own", dimless, 1.0)
    );
    
    surfaceScalarField phiNei
    (
        "phiNei",
        fvc::interpolate(U_, nei, "reconstruct(U)") & U_.mesh().Sf()
    );
    
    surfaceScalarField phiOwn
    (
        "phiOwn",
        fvc::interpolate(U_, own, "reconstruct(U)") & U_.mesh().Sf()
    );
    
    // Update interpolated nodes
    quadrature_.interpolateNodes();
    
    // Updated reconstructed moments
    quadrature_.momentsNei().update();
    quadrature_.momentsOwn().update();

    forAll(quadrature_.moments(), mI)
    {
        volUnivariateMoment& m = quadrature_.moments()[mI];
        
        // Old moment value
        const scalarField& psi0 = m.oldTime();

        // New moment value
        scalarField& psiIf = m;
        
        // Resetting moment
        psiIf = scalar(0);
        
        dimensionedScalar zeroPhi("zero", phiNei.dimensions(), 0.0);
         
        surfaceScalarField mFlux
        (
            quadrature_.momentsNei()[mI]*min(phiNei, zeroPhi) 
          + quadrature_.momentsOwn()[mI]*max(phiOwn, zeroPhi)
        );
         
         fvc::surfaceIntegrate(psiIf, mFlux);

         // Updating moment - First-order time integration here 
         psiIf = psi0 - psiIf*U_.time().deltaTValue();
    }
    
    quadrature_.updateQuadrature();
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceModels::univariatePopulationBalance::calcAggregation
(
    label order
)
{
    tmp<volScalarField> aSource
    (
        new volScalarField
        (
            IOobject
            (
                "aSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    volScalarField& aggregationSource = aSource();
       
    forAll(quadrature_.nodes(), pNode1I)
    {
        const volScalarNode& node1 = quadrature_.nodes()[pNode1I];
        
        const volScalarField& pWeight1 = node1.primaryWeight();
        
        forAll(node1.secondaryWeights(), sNode1I)
        {              
            
            const volScalarField& sWeight1 = node1.secondaryWeights()[sNode1I];
            
            const volScalarField& sAbscissa1 
                = node1.secondaryAbscissae()[sNode1I];
            
            forAll(quadrature_.nodes(), pNode2I)
            {
                const volScalarNode& node2 = quadrature_.nodes()[pNode2I];
                
                const volScalarField& pWeight2 = node2.primaryWeight();
                
                forAll(node2.secondaryWeights(), sNode2I)
                {
                    const volScalarField& sWeight2 
                        = node2.secondaryWeights()[sNode2I];

                    const volScalarField& sAbscissa2 
                        = node2.secondaryAbscissae()[sNode2I];
                    
                    tmp<volScalarField> aggInnerSum =
                        pWeight1*sWeight1*
                        (
                            pWeight2*sWeight2*
                            (
                                0.5*pow // Birth 
                                (
                                    pow3(sAbscissa1) + pow3(sAbscissa2), 
                                    order/3.0
                                )
                              - pow(sAbscissa1, order)
                            )*aggregationKernel_->Ka(sAbscissa1, sAbscissa2)
                        );
                                                
                    aggregationSource.dimensions().reset
                    (
                        aggInnerSum().dimensions()
                    );
                    
                    aggregationSource == aggregationSource + aggInnerSum();
                }
            }
        }
    }
    
    return aSource;
}

Foam::tmp<Foam::volScalarField> 
Foam::populationBalanceModels::univariatePopulationBalance::calcBreakup
(
    label order
)
{
    tmp<volScalarField> bSource
    (
        new volScalarField
        (
            IOobject
            (
                "bSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
        
    volScalarField& breakupSource = bSource();  
       
    forAll(quadrature_.nodes(), pNodeI)
    {
        const volScalarNode& node = quadrature_.nodes()[pNodeI];
        
        forAll(node.secondaryWeights(), sNodeI)
        {
            tmp<volScalarField> bSrc = node.primaryWeight()
                *node.secondaryWeights()[sNodeI]
                *breakupKernel_->Kb(node.secondaryAbscissae()[sNodeI])
                *(   
                    daughterDistribution_->mD                      //Birth
                    (
                        order, 
                        node.secondaryAbscissae()[sNodeI]
                    ) 
                  - pow(node.secondaryAbscissae()[sNodeI], order)   //Death
                 );
             
            breakupSource.dimensions().reset(bSrc().dimensions());
            breakupSource == breakupSource + bSrc;
        }
    }
    
    return bSource;
}

void Foam::populationBalanceModels::univariatePopulationBalance::solve()
{
    Info<< "Solving population balance equation.\n" << endl;

    // Advect moments with kinetic fluxes
    advectMoments();
    
    // Integrate source and diffusion terms
    forAll(quadrature_.moments(), mI)
    {
        volUnivariateMoment& m = quadrature_.moments()[mI];
        
        fvScalarMatrix moment
        (
            fvm::ddt(m) - fvc::ddt(m) 
        );
        
        if (aggregation_)
        {
            moment -= calcAggregation(m.order());
        }
        
        if (breakup_)
        {
            moment -= calcBreakup(m.order());
        }
        
        moment.relax();
        moment.solve();
    }
    
    quadrature_.updateQuadrature();
}


// ************************************************************************* //
