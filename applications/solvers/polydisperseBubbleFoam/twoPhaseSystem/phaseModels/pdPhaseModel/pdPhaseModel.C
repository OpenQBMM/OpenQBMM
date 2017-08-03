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
#include "twoPhaseSystem.H"
#include "fixedValueFvPatchFields.H"
#include "cyclicFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "emptyFvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "momentFieldSets.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pdPhaseModel::updateVelocity()
{
    // Correct mean velocity using the new velocity moments
    U_ =
        quadrature_.velocityMoments()[1]
       /Foam::max
        (
            quadrature_.moments()[1],
            residualAlpha_*rho()
        );

    U_.correctBoundaryConditions();
    phiPtr_() = fvc::flux(U_);
    alphaPhi_ = phiPtr_()*fvc::interpolate(*this);
    correctInflowOutflow(alphaPhi_);
    alphaRhoPhi_ = alphaPhi_*fvc::interpolate(rho());
}


Foam::tmp<Foam::volScalarField> Foam::pdPhaseModel::coalesenceSource
(
    const label& momentOrder
)
{
    tmp<volScalarField> cSource
    (
        new volScalarField
        (
            IOobject
            (
                "cSource",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar
            (
                "0",
                quadrature_.moments()[momentOrder].dimensions()/dimTime,
                0.0
            )
        )
    );

    if (!coalesence_)
    {
        return cSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, pNode1i)
    {
        const volScalarNode& node1 = nodes[pNode1i];
        const volScalarField& pWeight1 = node1.primaryWeight();
        const volScalarField& pAbscissa1 = node1.primaryAbscissa();

        forAll(nodes, pNode2i)
        {
            const volScalarNode& node2 = nodes[pNode2i];
            const volScalarField& pWeight2 = node2.primaryWeight();
            const volScalarField& pAbscissa2 = node2.primaryAbscissa();

            //- Diameter is used to calculate the coalesence kernel in place
            //  of the abscissa
            cSource.ref() +=
                pWeight1*
                (
                    pWeight2*
                    (
                        0.5*pow // Birth
                        (
                            pow3(pAbscissa1)
                          + pow3(pAbscissa2),
                            momentOrder/3.0
                        )
                      - pow(pAbscissa1, momentOrder)
                    )*fluid_.coalesence().Ka
                    (
                        ds_[pNode1i], ds_[pNode2i]
                    )
                );
        }
    }
    return cSource;
}


Foam::tmp<Foam::volScalarField> Foam::pdPhaseModel::breakupSource
(
    const label& momentOrder
)
{
    tmp<volScalarField> bSource
    (
        new volScalarField
        (
            IOobject
            (
                "bSource",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar
            (
                "0",
                quadrature_.moments()[momentOrder].dimensions()/dimTime,
                0.0
            )
        )
    );

    if (!breakup_)
    {
        return bSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(bSource(), celli)
    {
        forAll(nodes, pNodei)
        {
            const volScalarNode& node = nodes[pNodei];

            //- Diameter is used to calculate the breakup kernel in place
            //  of the abscissa
            bSource.ref()[celli] = bSource()[celli]
              + node.primaryWeight()[celli]
               *fluid_.breakup().Kb(ds_[pNodei][celli], celli)
               *(
                    fluid_.daughterDistribution().mD              //Birth
                    (
                        momentOrder,
                        node.primaryAbscissa()[celli]
                    )
                  - pow(node.primaryAbscissa()[celli], momentOrder)   //Death
                );
        }
    }

    return bSource;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdPhaseModel::pdPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid,phaseProperties,phaseName),
    pbeDict_
    (
        fluid.mesh().lookupObject<IOdictionary>("populationBalanceProperties")
    ),
    coalesence_(pbeDict_.lookup("coalesence")),
    breakup_(pbeDict_.lookup("breakup")),
    quadrature_(phaseName, fluid.mesh(), "RPlus"),
    nNodes_(quadrature_.nodes().size()),
    nMoments_(quadrature_.nMoments()),
    alphas_(nNodes_),
    Us_(quadrature_.velocities()),
    Vs_(nNodes_),
    ds_(nNodes_),
    maxD_("maxD", dimLength, phaseDict_),
    minD_("minD", dimLength, phaseDict_)
{
    if (nNodes_ == 1)
    {
        FatalErrorInFunction
            << "Polydisperse phase model selected, but only one node " << nl
            << "is used. Please use monodispersePhaseModel instead." << endl
            << exit(FatalError);
    }
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
                            phaseModel::name_,
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh(),
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
                            phaseModel::name_,
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedVector("zeroV", dimVelocity, Zero)
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
                            phaseModel::name_,
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh(),
                minD_
            )
        );
    }

    // Set alpha value based on moments
    volScalarField(*this) == quadrature_.moments()[1]/rho();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdPhaseModel::~pdPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pdPhaseModel::correct()
{
    quadrature_.updateAllQuadrature();

    d_ = dimensionedScalar("zero", dimLength, 0.0);

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
                    node.primaryWeight()*node.primaryAbscissa()/rho()
                   *(*this)/Foam::max
                    (
                        quadrature_.moments()[1]/rho(),
                        residualAlpha_
                    ),
                    dimensionedScalar("zero", dimless, 0.0)
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
                            dimensionedScalar("zero", dimMass, 0.0)
                        )*6.0
                       /(rho()*Foam::constant::mathematical::pi),
                        1.0/3.0
                    ),
                    minD_
                ),
                maxD_
            );

        d_ += alphas_[nodei]*ds_[nodei];
    }

    d_ /= Foam::max((*this), residualAlpha_);
    d_.max(minD_);
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
            dimensionedScalar zeroPhi("zero", phiPtr_().dimensions(), 0.0);

            // Create total flux field so that the individual fluxes can be
            // summed together
            volScalarField relativeDivVp
            (
                IOobject
                (
                    "relativeDivVp",
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fluid_.mesh(),
                dimensionedScalar("zero", m.dimensions()/dimTime, 0.0)
            );

            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
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
            }

            // Solve relative size moment transport equation
            fvScalarMatrix mEqn
            (
                fvm::ddt(m)
              + relativeDivVp
            );

            mEqn.relax();
            mEqn.solve();
        }

        forAll(quadrature_.velocityMoments(), mEqni)
        {
            volVectorField& Up = quadrature_.velocityMoments()[mEqni];
            dimensionedScalar zeroPhi("zero", phiPtr_().dimensions(), 0.0);

            // Create total flux field so that the individual fluxes can be
            // summed together
            volVectorField relativeDivPp
            (
                IOobject
                (
                    "relativeDivPp",
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fluid_.mesh(),
                dimensionedVector("zero", Up.dimensions()/dimTime, Zero)
            );

            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
                surfaceScalarField phiv("phiv", fvc::flux(Vs_[nodei]));

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
        }

        quadrature_.updateAllQuadrature();
        this->updateVelocity();
    }
    correct();
}

void Foam::pdPhaseModel::averageTransport(const PtrList<fvVectorMatrix>& AEqns)
{
    Info<< "Transporting moments with average velocity" << endl;

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
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
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
          + breakupSource(mEqni)
          + coalesenceSource(mEqni)
        );
        mEqn.relax();
        mEqn.solve();
    }

    forAll(quadrature_.velocityMoments(), mEqni)
    {
        dimensionedScalar zeroPhi("zero", phiPtr_().dimensions(), 0.0);
        volVectorField& Up = quadrature_.velocityMoments()[mEqni];

        volVectorField meanDivUbUp
        (
            IOobject
            (
                "meanDivUbUp",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
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
    correct();

    Info << "Solving for velocity abscissae" << endl;

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

        volScalarField alphaRhoi
        (
            "alphaRhoi",
            quadrature_.nodes()[nodei].primaryAbscissa()
           *quadrature_.nodes()[nodei].primaryWeight()
        );

        // Solve for velocities using acceleration terms
        fvVectorMatrix UsEqn
        (
            fvm::ddt(alphaRhoi, Us_[nodei])
          - alphaRhoi*fvc::ddt(Us_[nodei])
          + fvm::Sp(tauC*alphaRhoi, Us_[nodei])

         ==
            AEqns[nodei]
          + tauC*alphaRhoi*U_
        );

        UsEqn.relax();
        UsEqn.solve();
    }

    quadrature_.updateAllMoments();
    this->updateVelocity();

    // Update deviation velocity
    forAll(Vs_, nodei)
    {
        Vs_[nodei] = Us_[nodei] - U_;
    }
}

// ************************************************************************* //
