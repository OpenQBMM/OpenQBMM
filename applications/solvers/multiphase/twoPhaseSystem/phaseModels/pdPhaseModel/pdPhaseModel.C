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
    alphaPhi_ = fvc::interpolate(quadrature_.moments()[1]/rho())*phiPtr_();
    alphaRhoPhi_ = alphaPhi_*fvc::interpolate(rho());
}


Foam::tmp<Foam::volScalarField> Foam::pdPhaseModel::coalescenceSource
(
    const label momentOrder
)
{
    tmp<volScalarField> tmpCSource
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

    if (!coalescence_)
    {
        return tmpCSource;
    }

    volScalarField& cSource = tmpCSource.ref();
    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node1 = nodes[nodei];
        const volScalarField& weight1 = node1.primaryWeight();
        const volScalarField& abscissa1 = node1.primaryAbscissa();

        forAll(nodes, nodej)
        {
            const volScalarNode& node2 = nodes[nodej];
            const volScalarField& weight2 = node2.primaryWeight();
            const volScalarField& abscissa2 = node2.primaryAbscissa();

            //- Diameter is used to calculate the coalesence kernel in place
            //  of the abscissa, handled inside kernel
            cSource +=
                weight1
               *(
                    weight2
                   *(
                        0.5*pow // Birth
                        (
                            pow3(abscissa1)
                          + pow3(abscissa2),
                            momentOrder/3.0
                        )
                      - pow(abscissa1, momentOrder)
                    )*coalescenceKernel_.Ka(nodei, nodej)
                );
        }
    }
    return tmpCSource*pos(0.8 - *this);
}


Foam::tmp<Foam::volScalarField> Foam::pdPhaseModel::breakupSource
(
    const label momentOrder
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

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];

        //- Diameter is used to calculate the breakup kernel in place
        //  of the abscissa
        bSource.ref() +=
            node.primaryWeight()
           *breakupKernel_->Kb(nodei)
           *(
                daughterDistribution                       //Birth
                (
                    momentOrder,
                    node.primaryAbscissa()
                )
              - pow(node.primaryAbscissa(), momentOrder)   //Death
            );
    }

    return bSource*pos(0.8 - *this);
}


Foam::tmp<Foam::volScalarField> Foam::pdPhaseModel::daughterDistribution
(
    const label momentOrder,
    const volScalarField& abscissa
)
{
    tmp<volScalarField> tmpDaughterDist
    (
        new volScalarField
        (
            IOobject
            (
                "daughterDist",
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
                pow(dimMass, momentOrder),
                0.0
            )
        )
    );
    volScalarField& daughterDist = tmpDaughterDist.ref();

    forAll(daughterDist, celli)
    {
        daughterDist[celli] = daughterDistribution_->mD
        (
            momentOrder,
            abscissa[celli]
        );
    }

    return tmpDaughterDist;
}


void Foam::pdPhaseModel::solveSourceOde()
{
    if (!ode_)
    {
        forAll(quadrature_.moments(), mI)
        {
            quadrature_.moments()[mI] +=
                U_.mesh().time().deltaT()
               *(
                    coalescenceSource(mI) + breakupSource(mI)
                );
        }
        quadrature_.updateQuadrature();
        quadrature_.updateAllMoments();
        return;
    }

    PtrList<volScalarField> k1(nMoments_);
    PtrList<volScalarField> k2(nMoments_);
    PtrList<volScalarField> k3(nMoments_);
    PtrList<volScalarField> momentsOld(nMoments_);

    forAll(momentsOld, mI)
    {
        k1.set
        (
            mI,
            new volScalarField
            (
                IOobject
                (
                    "k1",
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar
                (
                    "k1",
                    quadrature_.moments()[mI].dimensions(),
                    0.0
                )
            )
        );

        k2.set
        (
            mI,
            new volScalarField
            (
                IOobject
                (
                    "k2",
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar
                (
                    "k2",
                    quadrature_.moments()[mI].dimensions(),
                    0.0
                )
            )
        );

        k3.set
        (
            mI,
            new volScalarField
            (
                IOobject
                (
                    "k3",
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar
                (
                    "k3",
                    quadrature_.moments()[mI].dimensions(),
                    0.0
                )
            )
        );

        momentsOld.set
        (
            mI,
            new volScalarField(quadrature_.moments()[mI])
        );
    }
    quadrature_.updateQuadrature();

    // Read current deltaT
    dimensionedScalar dt0 = U_.mesh().time().deltaT();


    // Calculate k1 for all moments
    forAll(momentsOld, mI)
    {
        k1[mI] = dt0*(coalescenceSource(mI) + breakupSource(mI));
        quadrature_.moments()[mI] == momentsOld[mI] + k1[mI];
    }
    quadrature_.updateQuadrature();

    // Calculate k2 for all moments
    forAll(momentsOld, mI)
    {
        k2[mI] = dt0*(coalescenceSource(mI) + breakupSource(mI));
        quadrature_.moments()[mI] == momentsOld[mI] + (k1[mI] + k2[mI])/4.0;
    }
    quadrature_.updateQuadrature();

    // calculate k3 and new moments for all moments
    forAll(momentsOld, mI)
    {
        k3[mI] = dt0*(coalescenceSource(mI) + breakupSource(mI));

        // Second order accurate, k3 only used for error estimation
        quadrature_.moments()[mI] ==
            momentsOld[mI]
          + (k1[mI] + k2[mI] + 4.0*k3[mI])/6.0;
    }
    quadrature_.updateQuadrature();

    // Because velocity moments only change due to change in size abscissae
    // from break up and coalescence, the velocity moments are simply updated
    // to include this size change. This eliminated calculating source terms
    // twice.
    quadrature_.updateAllMoments();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdPhaseModel::pdPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    pbeDict_
    (
        IOobject
        (
            "populationBalanceProperties",
            fluid.mesh().time().constant(),
            fluid.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    ode_(pbeDict_.lookupOrDefault("ode", false)),
    coalescence_(pbeDict_.lookup("coalescence")),
    breakup_(pbeDict_.lookup("breakup")),
    quadrature_(phaseName, fluid.mesh(), "RPlus"),
    ddtM1_
    (
        IOobject
        (
            IOobject::groupName("ddtM1", phaseName),
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("0", dimDensity/dimTime, 0.0)
    ),
    nNodes_(quadrature_.nodes().size()),
    nMoments_(quadrature_.nMoments()),
    alphas_(nNodes_),
    Us_(quadrature_.velocities()),
    Vs_(nNodes_),
    ds_(nNodes_),
    maxD_("maxD", dimLength, phaseDict_),
    minD_("minD", dimLength, phaseDict_),
    coalescenceKernel_
    (
        pbeDict_.subDict("coalescenceKernel"),
        fluid.mesh()
    ),
    breakupKernel_
    (
        Foam::bubbleBreakupKernel::New
        (
            pbeDict_.subDict("breakupKernel"),
            fluid.mesh()
        )
    ),
    daughterDistribution_
    (
        Foam::populationBalanceSubModels::daughterDistribution::New
        (
            pbeDict_.subDict("daughterDistribution")
        )
    )
{
    this->d_.writeOpt() = IOobject::AUTO_WRITE;

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
    d_ = dimensionedScalar("zero", dimLength, 0.0);

    volScalarField scale
    (
        (*this)
       /Foam::max
        (
            quadrature_.moments()[1]/rho(),
            residualAlpha_
        )
    );

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
                node.primaryWeight()*node.primaryAbscissa()/rho()*scale;
            alphas_[nodei].max(0);
            alphas_[nodei].min(1);
        }

        //  Calculate bubble diameter based on bubble mass (abscissa)
        ds_[nodei] =
            Foam::min
            (
                Foam::max
                (
                    Foam::pow
                    (
                        node.primaryAbscissa()*6.0
                       /(rho()*Foam::constant::mathematical::pi)
                      + dimensionedScalar("smallVolume", dimVolume, SMALL),
                        1.0/3.0
                    ),
                    minD_
                ),
                maxD_
            );

        if (nNodes_ > 1)
        {
            d_ += alphas_[nodei]*ds_[nodei];
        }
    }

    if (nNodes_ > 1)
    {
        d_ /= Foam::max((*this), residualAlpha_);
    }
    else
    {
        d_ = ds_[0];
    }
    d_.max(minD_);
}


void Foam::pdPhaseModel::relativeTransport()
{
    // Do nothing if only mean is used
    if (nNodes_ == 1)
    {
        return;
    }

    Info<< "Transporting moments based on relative flux" << endl;

    quadrature_.interpolateNodes();
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_.nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_.nodesNei();

    // Transport moments with relative flux
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
    ddtM1_ = fvc::ddt(quadrature_.moments()[1]);

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
    correct();

    quadrature_.interpolateNodes();
}


void Foam::pdPhaseModel::averageTransport(const PtrList<fvVectorMatrix>& AEqns)
{
    // Update moments based source terms for breakup and coalescence
    solveSourceOde();

    // Mean moment advection
    Info<< "Transporting moments with average velocity" << endl;
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_.nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_.nodesNei();

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
        );
        mEqn.relax();
        mEqn.solve();
    }

    if (nNodes_ == 1)
    {
        forAll(quadrature_.velocityMoments(), mi)
        {
            quadrature_.velocityMoments()[mi] = U_*quadrature_.moments()[mi];
            quadrature_.velocityMoments()[mi].correctBoundaryConditions();
        }

        quadrature_.updateAllQuadrature();
        return;
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

    // Solve for velocity abscissa directly since the momentum exchange
    // terms do not change the mass
    Info << "Solving for velocity abscissae" << endl;

    forAll(Us_, nodei)
    {
        //  Colisional time, forces velocities towards mean in the case of
        //  high volume fractions
        //  Could be replaced by radial distribution function
        volScalarField tauC
        (
            "tauC",
            (0.5 + 0.5*tanh(((*this) - 0.63)/0.01))*HUGE
        );

        tauC.dimensions().reset(inv(dimTime));

        volScalarField alphaRhoi
        (
            "alphaRhoi",
            alphas_[nodei]*rho()
        );

        // Solve for velocities using acceleration terms
        fvVectorMatrix UsEqn
        (
            alphaRhoi*fvm::ddt(Us_[nodei])
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
    correct();

    // Update deviation velocity
    forAll(Vs_, nodei)
    {
        Vs_[nodei] = Us_[nodei] - U_;
    }
}

// ************************************************************************* //
