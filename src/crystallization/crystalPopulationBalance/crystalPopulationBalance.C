/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
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

#include "crystalPopulationBalance.H"
#include "rhoReactionThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(crystalPopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        crystalPopulationBalance,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::crystalPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    univariatePDFTransportModel(name, dict, phi.mesh(), phi, "RPlus"),
    populationBalanceModel(name, dict, phi),
    odeType(phi.mesh(), dict),
    fixedSize_(dict.lookupOrDefault("fixedSize", false)),
    minM0_(readScalar(dict.lookup("minM0"))),
    nucleation_(dict.lookupOrDefault("nucleation", false)),
    growth_(dict.lookupOrDefault("growth", false)),
    aggregation_(dict.lookupOrDefault("aggregation", false)),
    breakup_(dict.lookupOrDefault("breakup", false)),
    deposition_(dict.lookupOrDefault("deposition", false)),
    saturation_(dict.lookup("saturation")),
    speciesCoupled_(dict.lookupOrDefault("speciesCoupled", false)),
    nucleationModel_(),
    growthModel_(),
    diffusionModel_
    (
        Foam::populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    aggregationKernel_(),
    breakupKernel_(),
    depositionModel_(),
    saturationModel_(),
    rhop_
    (
        dict.lookupOrDefault
        (
            "rhop",
            dimensionedScalar("rhop", dimDensity, 1000.0)
        )
    ),
    shapeFactor_(
        dict.lookupOrDefault
        (
            "shapeFactor",
            dimensionedScalar("shapeFactor", dimless, constant::mathematical::pi/6.0)
        )
    ),
    Pr_(
        dict.lookupOrDefault
        (
            "Pr",
            dimensionedScalar("Pr", dimless, 1.0)
        )
    ),
    Prt_(
        dict.lookupOrDefault
        (
            "Prt",
            dimensionedScalar("Prt", dimless, 1.0)
        )
    ),
    Sct_(
        dict.lookupOrDefault
        (
            "Sct",
            dimensionedScalar("Sct", dimless, 0.7)
        )
    ),
    molarDiff_(
        dict.lookupOrDefault
        (
            "molarDiff",
            dimensionedScalar("molarDiff", dimViscosity, 1.0e-5)
        )
    ),
    L10_
    (
        IOobject
        (
            "L10.crystal",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimLength,
        fvPatchFieldBase::zeroGradientType()
    ),
    L32_
    (
        IOobject
        (
            "L32.crystal",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimLength,
        fvPatchFieldBase::zeroGradientType()
    ),
    d_
    (
        IOobject
        (
            "d.crystal",
            phi.mesh().time().timeName(),
            phi.mesh(),                    
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi.mesh(),
        dimLength,
        fvPatchFieldBase::zeroGradientType()
    ),
    SYact_
    (
        IOobject
        (
            "SYact",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    d_nucleation_
    (
        "d_nucleation",
        dimLength,
        dict.lookupOrDefault<scalar>("d_nucleation", dict.lookupOrDefault<scalar>("deltaWidth", 1.0e-6))
    )
{
    if (quadrature_.nodes().size() && quadrature_.nodes()[0].extended())
    {
        FatalErrorInFunction
            << "crystalPopulationBalance supports non-extended quadrature only."
            << abort(FatalError);
    }

    if (saturation_ )
    {
        saturationModel_.reset
        (
            solutionSaturationModel::New
            (
                dict.subDict("saturationModel"),
                phi_.mesh()
            )            
        );
    }

    if (nucleation_)
    {
        nucleationModel_ =
            Foam::populationBalanceSubModels::nucleationModel::New
            (
                dict.subDict("nucleationModel"),
                phi_.mesh()
            );
    }

    if (growth_)
    {
        growthModel_ =
            Foam::populationBalanceSubModels::growthModel::New
            (
                dict.subDict("growthModel"),
                phi_.mesh()
            );
    }


    if (aggregation_)
    {
        aggregationKernel_ =
            Foam::populationBalanceSubModels::aggregationKernel::New
            (
                dict.subDict("aggregationKernel"),
                phi_.mesh()
            );
    }

    if (breakup_)
    {
        breakupKernel_ =
            Foam::populationBalanceSubModels::breakupKernel::New
            (
                dict.subDict("breakupKernel"),
                phi_.mesh()
            );
    }

    if (deposition_)
    {
        depositionModel_ = 
             Foam:: populationBalanceSubModels::depositionModel::New
            (
                dict.subDict("depositionModel"),
                phi_.mesh()
            );
    }

    updateCharacteristicLengths();

    if (speciesCoupled() && !growthModel_)
    {
        FatalErrorInFunction
            << "speciesCoupled requires an active growthModel."
            << abort(FatalError);
    }

    if (fixedSize_)
    {
        const    volScalarMomentFieldSet& moments = quadrature_.moments();
        PtrList<volScalarNode>& nodes = quadrature_.nodes();
        PtrList<volScalarField>& absc = nodes[0].primaryAbscissae();
        volScalarField& weig = nodes[0].primaryWeight();
        absc[0] = d_nucleation_.value();
        weig = moments[0];
        L10_ = d_nucleation_;
        L32_ = d_nucleation_;
        d_ = d_nucleation_;
        L10_.correctBoundaryConditions();
        L32_.correctBoundaryConditions();
        d_.correctBoundaryConditions();
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::~crystalPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



void
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::subModelsPreUpdate()
{
    if (aggregation_)
    {
        aggregationKernel_->preUpdate();
    }

    if (breakup_)
    {
        breakupKernel_->preUpdate();
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::implicitMomentSource
(
    const volScalarMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}


void
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::updateCellMomentSource(const label)
{}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
)
{

    scalar source(0);

    if (nucleation_)
    {
        source += 
            nucleationModel_->nucleationSource
            (
                momentOrder[0], 
                celli
            );
    }
    
    if (growth_)
    {
        source +=
            growthModel_->phaseSpaceConvection
            (
                momentOrder,
                celli,
                quadrature
            );
    }

    
    if (aggregation_)
    {
        source +=
            aggregationKernel_->aggregationSource
            (
                momentOrder,
                celli,
                quadrature,
                environment
            );
    }

    if (breakup_)
    {
        source +=
            breakupKernel_->breakupSource
            (
                momentOrder,
                celli,
                quadrature
            );
    }


//    if(celli == 0) 
//        Info << "source celli0 momentOrder: " << momentOrder[0] << tab <<  source << endl ;

    return source;
    
}



Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::realizableCo() const
{
    return univariatePDFTransportModel::realizableCo();
}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::CoNum() const
{
    return 0.0;
}


bool
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::solveMomentSources() const
{
    if (nucleation_ ||  growth_ || aggregation_ ||   deposition_)
    {
        return odeType::solveSources_;
    }

    return false;
}


bool
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::solveMomentOde() const
{
    return odeType::solveOde_;
}


void
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::checkCorrectMoments()
{

    volScalarMomentFieldSet& moments = quadrature_.moments();

    const scalar cutOffSize = d_nucleation_.value();
    const scalar numTol = 1;

    forAll(moments[0], celli)
    {
        scalar m0 = moments[0][celli];
        if (m0 <= numTol)
        {
            forAll(moments, mi)
            {
                moments[mi][celli] = 0.0;
            }
        }
        else 
        {
            scalar m1 = moments[1][celli];
            scalar meanL = m1 / (m0 + SMALL);

            if (meanL < cutOffSize)
            {
                forAll(moments, mi)
                {
                    moments[mi][celli] = 0.0;
                }
            }
        }
    }
    forAll(moments, mi)
        moments[mi].correctBoundaryConditions();
}


void
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::updateCharacteristicLengths()
{
    if (fixedSize_)
    {
        return;
    }

    const volScalarMomentFieldSet& moments = quadrature_.moments();

    forAll(L10_, celli)
    {
        const scalar m0 = moments[0][celli];
        const scalar m1 = moments.size() > 1 ? moments[1][celli] : 0.0;
        const scalar m2 = moments.size() > 2 ? moments[2][celli] : 0.0;
        const scalar m3 = moments.size() > 3 ? moments[3][celli] : 0.0;

        L10_[celli] = m1/(m0 + SMALL);
        L32_[celli] =
            moments.size() > 3
          ? m3/(m2 + SMALL)
          : L10_[celli];

        d_[celli] = L32_[celli];
    }

    L10_.correctBoundaryConditions();
    L32_.correctBoundaryConditions();
    d_.correctBoundaryConditions();
}

void
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::calcSpeciesTransfer()
{
    SYact_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);


    if(speciesCoupled())
    {

        labelList momentOrder(1);
        momentOrder[0] = 3;
        
        // phaseSpaceConvection returns dimensionless dm3/dt
        
        scalar srcFactor = shapeFactor_.value() * rhop_.value();
        
        static bool firstTime = true;
        if (firstTime)
        {
            Info<< "calcSpeciesTransfer parameters:" << nl
                << "  shapeFactor = " << shapeFactor_.value() << nl
                << "  rhop = " << rhop_.value() << " kg/m³" << nl
                << "  srcFactor = " << srcFactor << " kg" << endl;
            firstTime = false;
        }

        scalar maxSource = 0.0;
        scalar maxSYact = 0.0;
        const volScalarMomentFieldSet& moments = quadrature_.moments();
        const volScalarField& C = phi_.mesh().lookupObject<volScalarField>("C");
        const scalar dt = phi_.mesh().time().deltaTValue();
        const scalar invDt = 1.0 / (dt + SMALL);
        const scalar maxSYact_abs = 100.0; 
        

        
        forAll(d_, celli)
        {
            scalar source = growthModel_->phaseSpaceConvection
            (
                momentOrder,
                celli,
                quadrature_
            );

            if(nucleation_)
                source += nucleationModel_->nucleationSource(momentOrder[0], celli);

            scalar SYact_final = 0.0;

            // calculate the source term
            scalar SYact_calc = srcFactor * source;

            // dissolution limit logic (for SYact_calc < 0)
            if (SYact_calc < -SMALL)
            {
                // check A: whether there are particles (using M0 to judge)
                // if there are few particles, do not allow dissolution
                if (moments[0][celli] < 1)
                {
                    SYact_calc = 0.0;
                }
                else
                {
                    // check B: mass conservation limit
                    // existing mass [kg/m3] = srcFactor * M3
                    scalar existingMass = srcFactor * moments[3][celli];
                    
                    scalar maxDissolutionRate = existingMass * invDt;
                    
                    if (mag(SYact_calc) > maxDissolutionRate)
                    {
                        SYact_final = -0.9 * maxDissolutionRate;
                    }
                    else
                    {
                        SYact_final = SYact_calc;
                    }
                } 
            }
            else
            {
                const scalar maxSYact_local = 0.9 * max(C[celli], SMALL) * invDt;

                SYact_final = min(SYact_calc, min(maxSYact_local, maxSYact_abs));
            }   

            SYact_[celli] = SYact_final;

            maxSource = max(maxSource, mag(source));
            maxSYact = max(maxSYact, mag(SYact_[celli]));
            
            // warning logic (optional, suggest only alarm when growth is limited, dissolution is usually normal)
            if (mag(SYact_calc) > mag(SYact_final) + SMALL && phi_.mesh().time().timeIndex() % 100 == 0)
            {
                Info<< "  max |source| (dm3/dt) = " << maxSource << " 1/s" << nl
                    << "  max |SYact| = " << maxSYact << " kg/(m³·s)" << endl;
            }
        
        }

    }
    SYact_.correctBoundaryConditions();

}


void
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::explicitMomentSource()
{
    odeType::solve(quadrature_, 0);

    if (deposition_)
    {
        depositionModel_->correct(quadrature_);
    }

}


void
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::solve()
{

    this->subModelsPreUpdate();

    // moment advection and diffusion  
    this->momentAdvection_().update();

    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());

    // Solve moment transport equations
    forAll(quadrature_.moments(), mEqni)
    {
        
        if(fixedSize_ && mEqni > 0) continue;

        volScalarMoment& m = quadrature_.moments()[mEqni];

        fvScalarMatrix mEqn
        (
            fvm::ddt(m)
            + momentAdvection_().divMoments()[mEqni]
            ==
            implicitMomentSource(m)
        );

        if(! (diffusionModel_->type() == "none"))
            mEqn.relax();

        mEqn.solve();

        m.correctBoundaryConditions();
    }

    if(! fixedSize_)
    {

        this->checkCorrectMoments();

        quadrature_.updateQuadrature();

        if (solveMomentSources())
        {
            this->explicitMomentSource();
        }

        updateCharacteristicLengths();

        calcSpeciesTransfer();

    }
    else
    {

        PtrList<volScalarNode>& nodes = quadrature_.nodes();
        nodes[0].primaryWeight() = quadrature_.moments()[0];
    
    }


}


bool 
Foam::PDFTransportModels::populationBalanceModels::crystalPopulationBalance
::readIfModified()
{
    odeType::read
    (
        populationBalanceProperties_
    );
    
    return true;
}


// ************************************************************************* //
