/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "pdPhaseModel.H"
#include "fvMatrix.H"

#include "fixedValueFvPatchFields.H"
#include "cyclicFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "emptyFvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdPhaseModel::pdPhaseModel
(
    const fvMesh& mesh,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(mesh,phaseProperties,phaseName),
    quadrature_(phaseName, mesh_, "RPlus"),
    nNodes_(quadrature_.nodes().size()),
    nMoments_(quadrature_.nMoments()),
    alphas_(nNodes_),
    Us_(quadrature_.velocities()),
    Vs_(nNodes_),
    ds_(nNodes_),
    maxD_("maxD", dimLength, phaseDict_),
    minD_("minD", dimLength, phaseDict_)
{
    wordList phiTypes
    (
        U_.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    forAll(U_.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
        )
        {
            phiTypes[i] = fixedValueFvsPatchScalarField::typeName;
        }
    }

    forAll(alphas_, nodei)
    {
        alphas_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "alpha",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("alpha", dimless, 0.0)
            )
        );
        Vs_.set
        (
            nodei,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "V",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("V", dimVelocity, Zero)
            )
        );
        ds_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "d",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("d", dimLength, 0.0)
            )
        );
    }

    // Set alpha value based on moments
    *this == quadrature_.moments()[1]/rho_;

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdPhaseModel::~pdPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pdPhaseModel::correct()
{
    quadrature_.updateAllQuadrature();

    d_ = dimensionedScalar("zero",dimLength, 0.0);

    forAll(quadrature_.nodes(), nodei)
    {
        const volScalarNode& node = quadrature_.nodes()[nodei];

        // Set alpha values such that the moment.1 is equal to the bounded
        // alpha
        if (nNodes_ == 1)
        {
            alphas_[nodei] = *this;
        }
        else
        {
            alphas_[nodei] =
                Foam::max
                (
                    node.primaryWeight()*node.primaryAbscissa()/rho_
                   *(*this)/Foam::max
                    (
                        quadrature_.moments()[1]/rho_,
                        residualAlpha_
                    ),
                    dimensionedScalar("0",dimless,0.0)
                );
        }

        //  Calculate bubble diameter based on bubble mass (abscissa)
        ds_[nodei] =
            Foam::min
            (
                Foam::max
                (
                    Foam::pow
                    (
                        Foam::max
                        (
                            node.primaryAbscissa(),
                            dimensionedScalar("0", dimMass, 0.0)
                        )*6.0
                       /(rho_*Foam::constant::mathematical::pi),
                        1.0/3.0
                    ),
                    minD_
                ),
                maxD_
            );
        d_ += alphas_[nodei]*ds_[nodei];
    }
    d_ /= Foam::max((*this), residualAlpha_);
}


void Foam::pdPhaseModel::relativeTransport()
{
    Info<< "Transporting moments based on relative flux" << endl;

    quadrature_.interpolateNodes();
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_.nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_.nodesNei();

    // Transport moments with relative flux only if polydisperse
    if (nNodes_ > 1)
    {
        forAll(quadrature_.moments(), mEqni)
        {
            volScalarField& m = quadrature_.moments()[mEqni];
            volVectorField& Up = quadrature_.velocityMoments()[mEqni];
            dimensionedScalar zeroPhi("zero", phiPtr_().dimensions(), 0.0);

            // Create total flux field so that the individual fluxes can be summed
            // together
            volScalarField relativeDivVp
            (
                IOobject
                (
                    "relativeDivVp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedScalar("zero", m.dimensions()/dimTime, 0.0)
            );

            volVectorField relativeDivPp
            (
                IOobject
                (
                    "relativeDivPp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedVector("zero", Up.dimensions()/dimTime, Zero)
            );

            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
                Vs_[nodei] = Us_[nodei] - U_;
                surfaceScalarField phiv("phiv", fvc::flux(Vs_[nodei]));

                // Calculate size moment flux
                surfaceScalarField rFluxVp
                (
                    nodesNei[nodei].primaryWeight()
                   *(
                        pow
                        (
                            nodesNei[nodei].primaryAbscissa(),
                            mEqni
                        )
                    )*Foam::min(phiv, zeroPhi)
                  + nodesOwn[nodei].primaryWeight()
                   *pow
                    (
                        nodesOwn[nodei].primaryAbscissa(),
                        mEqni
                    )*Foam::max(phiv, zeroPhi)
                );

                relativeDivVp += fvc::surfaceIntegrate(rFluxVp);

                // Calculate velocity moment flux
                surfaceVectorField rFluxPp
                (
                    "rFluxPp",
                    quadrature_.velocitiesNei()[nodei]
                   *nodesNei[nodei].primaryWeight()
                    *(
                        pow
                        (
                            nodesNei[nodei].primaryAbscissa(),
                            mEqni
                        )
                    )*Foam::min(phiv, zeroPhi)
                  + quadrature_.velocitiesOwn()[nodei]
                   *nodesOwn[nodei].primaryWeight()
                   *pow
                    (
                        nodesOwn[nodei].primaryAbscissa(),
                        mEqni
                    )*Foam::max(phiv, zeroPhi)
                );

                relativeDivPp += fvc::surfaceIntegrate(rFluxPp);
            }


            // Solve relative velocity moment transport equation
            fvVectorMatrix UpEqn
            (
                fvm::ddt(Up)
              + relativeDivPp
            );

            UpEqn.relax();
            UpEqn.solve();


            // Solve relative size moment transport equation
            fvScalarMatrix mEqn
            (
                fvm::ddt(m)
              + relativeDivVp
            );

            mEqn.relax();
            mEqn.solve();
        }
        quadrature_.updateAllQuadrature();

        // Update mean velocity based on new velocity moments
        U_ =
            quadrature_.velocityMoments()[1]
           /Foam::max
            (
                quadrature_.moments()[1],
                residualAlpha_*rho_
            );

        U_.correctBoundaryConditions();
        phiPtr_() == fvc::flux(U_);
    }
}

void Foam::pdPhaseModel::averageTransport(const PtrList<fvVectorMatrix>& AEqns)
{
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_.nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_.nodesNei();

    quadrature_.interpolateNodes();

    forAll(quadrature_.moments(), mEqni)
    {
        volScalarField& m = quadrature_.moments()[mEqni];

        dimensionedScalar zeroPhi("zero", phiPtr_().dimensions(), 0.0);
        volScalarField meanDivUbMp
        (
            IOobject
            (
                "meanDivUbMp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", m.dimensions()/dimTime, Zero)
        );

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            // Update average size moment flux
            surfaceScalarField aFluxMp
            (
                "aFluxMp",
                nodesNei[nodei].primaryWeight()
               *(
                    pow
                    (
                        nodesNei[nodei].primaryAbscissa(),
                        mEqni
                    )
                )*Foam::min(phiPtr_(), zeroPhi)
              + nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissa(),
                    mEqni
                )*Foam::max(phiPtr_(), zeroPhi)
            );

            meanDivUbMp += fvc::surfaceIntegrate(aFluxMp);
        }

        // Solve average size moment transport
        fvScalarMatrix mEqn
        (
            fvm::ddt(m)
          - fvc::ddt(m)
          + meanDivUbMp
        );

        mEqn.relax();
        mEqn.solve();
    }
    correct();

    // If momodisperse, use mean velocity to construct velocity moments
//     if(nNodes_ == 1)
//     {
//         forAll(quadrature_.velocityMoments(), mi)
//         {
//             quadrature_.velocityMoments()[mi] = U_*quadrature_.moments()[mi];
//             quadrature_.velocityMoments()[mi].correctBoundaryConditions();
//         }
//
//         quadrature_.updateAllQuadrature();
//     }
//     else
    {
        forAll(quadrature_.velocityMoments(), mEqni)
        {
            dimensionedScalar zeroPhi("zero", phiPtr_().dimensions(), 0.0);
            volVectorField& Up = quadrature_.velocityMoments()[mEqni];

            volVectorField meanDivUbUp
            (
                IOobject
                (
                    "meanDivUbUp",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedVector("zero", Up.dimensions()/dimTime, Zero)
            );

            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
                // Update average velocity moment flux
                surfaceVectorField aFluxUp
                (
                    "aFluxUp",
                    quadrature_.velocitiesNei()[nodei]
                   *nodesNei[nodei].primaryWeight()
                    *(
                        pow
                        (
                            nodesNei[nodei].primaryAbscissa(),
                            mEqni
                        )
                    )*Foam::min(phiPtr_(), zeroPhi)
                  + quadrature_.velocitiesOwn()[nodei]
                   *nodesOwn[nodei].primaryWeight()
                   *pow
                    (
                        nodesOwn[nodei].primaryAbscissa(),
                        mEqni
                    )*Foam::max(phiPtr_(), zeroPhi)
                );

                meanDivUbUp += fvc::surfaceIntegrate(aFluxUp);
            }


            // Solve average velocity moment transport Equation
            fvVectorMatrix UpEqn
            (
                fvm::ddt(Up)
              - fvc::ddt(Up)
              + meanDivUbUp
            );

            UpEqn.relax();
            UpEqn.solve();
        }
        quadrature_.updateAllQuadrature();

        // Solve for velocity abscissa directly since the momentum exchange
        //  terms do not change the mass
        forAll(Us_, nodei)
        {
            //  Colisional time, forces velocities towards mean in the case of
            //  high volume fractions
            volScalarField tauC
            (
                "tauC",
                (0.5 + 0.5*tanh(((*this) - 0.63)/0.01))*HUGE
            );
            tauC.dimensions().reset(inv(dimTime));

            // Solve for velocities using acceleration terms
            fvVectorMatrix UsEqn
            (
                fvm::ddt(Us_[nodei])
              - fvc::ddt(Us_[nodei])
              + fvm::Sp(tauC, Us_[nodei])

            ==
                AEqns[nodei]
              + tauC*U_
            );

            UsEqn.relax();
            UsEqn.solve();
        }
        updateMoments();
    }
}

void Foam::pdPhaseModel::updateMoments()
{
    quadrature_.updateAllMoments();

    // Correct mean velocity using the new velocity moments
    U_ =
        quadrature_.velocityMoments()[1]
       /Foam::max
        (
            quadrature_.moments()[1],
            residualAlpha_*rho_
        );

    U_.correctBoundaryConditions();
    phiPtr_() == fvc::flux(U_);
}



// ************************************************************************* //
