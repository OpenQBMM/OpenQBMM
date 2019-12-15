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

#include "polydispersePhaseModel.H"
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
#include "vectorList.H"
#include "fixedFaceFvPatchScalarField.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polydispersePhaseModel::updateVelocity()
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
    alphaPhi_ = fvc::interpolate(*this)*phiPtr_();
    alphaRhoPhi_ = fvc::interpolate(rho())*alphaPhi_;
}


Foam::scalar Foam::polydispersePhaseModel::coalescenceSource
(
    const label celli,
    const label momentOrder
)
{
    scalar cSource = 0.0;
    if (!coalescence_)
    {
        return cSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node1 = nodes[nodei];
        scalar weight1 = node1.primaryWeight()[celli];
        scalar abscissa1 = Foam::max(node1.primaryAbscissae()[0][celli], SMALL);
        scalar n1 = node1.n(celli, weight1, abscissa1);
        scalar d1 = node1.d(celli, abscissa1);

        forAll(nodes, nodej)
        {
            const volScalarNode& node2 = nodes[nodej];
            scalar weight2 = node2.primaryWeight()[celli];
            scalar abscissa2 = Foam::max(node2.primaryAbscissae()[0][celli], SMALL);

            scalar n2 = node2.n(celli, weight2, abscissa2);
            scalar d2 = node2.d(celli, abscissa2);
            vector Ur = Us_[nodei][celli] - Us_[nodej][celli];

            //- Diameter is used to calculate the coalesence kernel in place
            //  of the abscissa, handled inside kernel
            cSource +=
                0.5*n1
               *(
                    n2
                   *(
                        pow // Birth
                        (
                            (abscissa1) + (abscissa2),
                            momentOrder
                        )
                      - pow(abscissa1, momentOrder)
                      - pow(abscissa2, momentOrder)
                    )
                )*coalescenceKernel_->Ka(d1, d2, Ur, celli);
        }
    }
    return cSource;
}

Foam::vector Foam::polydispersePhaseModel::coalescenceSourceU
(
    const label celli,
    const label momentOrder
)
{
    vector cSource = Zero;
    if (!coalescence_ || momentOrder == 1)
    {
        return cSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node1 = nodes[nodei];
        scalar weight1 = node1.primaryWeight()[celli];
        scalar abscissa1 = Foam::max(node1.primaryAbscissae()[0][celli], SMALL);
        scalar n1 = node1.n(celli, weight1, abscissa1);
        scalar d1 = node1.d(celli, abscissa1);

        forAll(nodes, nodej)
        {
            const volScalarNode& node2 = nodes[nodej];
            scalar weight2 = node2.primaryWeight()[celli];
            scalar abscissa2 = Foam::max(node2.primaryAbscissae()[0][celli], SMALL);

            scalar n2 = node2.n(celli, weight2, abscissa2);
            scalar d2 = node2.d(celli, abscissa2);
            vector Ur = Us_[nodei][celli] - Us_[nodej][celli];

            //- Diameter is used to calculate the coalesence kernel in place
            //  of the abscissa, handled inside kernel
            cSource +=
                0.5*n1*Us_[nodei][celli]
               *(
                    n2
                   *(
                        pow // Birth
                        (
                            (abscissa1) + (abscissa2),
                            momentOrder
                        )
                      - pow(abscissa1, momentOrder)
                      - pow(abscissa2, momentOrder)
                    )*coalescenceKernel_->Ka(d1, d2, Ur, celli)
                );
        }
    }

    return cmptMultiply(cSource, validDirections_);
}


Foam::scalar Foam::polydispersePhaseModel::breakupSource
(
    const label celli,
    const label momentOrder
)
{
    scalar bSource = 0.0;
    if (!breakup_)
    {
        return bSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];
        scalar weight = node.primaryWeight()[celli];
        scalar abscissa = Foam::max(node.primaryAbscissae()[0][celli], SMALL);
        scalar d = node.d(celli, abscissa);
        scalar n = node.n(celli, weight, abscissa);

        //- Diameter is used to calculate the breakup kernel in place
        //  of the abscissa
        bSource +=
            n
           *breakupKernel_->Kb(d, celli)
           *breakupKernel_->massNodeSource  //Birth and death
            (
                abscissa,
                momentOrder
            );
    }
    return bSource;
}


Foam::vector Foam::polydispersePhaseModel::breakupSourceU
(
    const label celli,
    const label momentOrder
)
{
    vector bSource = Zero;
    if (!breakup_  || momentOrder == 1)
    {
        return bSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];
        scalar weight = node.primaryWeight()[celli];
        scalar abscissa = Foam::max(node.primaryAbscissae()[0][celli], SMALL);
        scalar d = node.d(celli, abscissa);
        scalar n = node.n(celli, weight, abscissa);

        //- Diameter is used to calculate the breakup kernel in place
        //  of the abscissa
        bSource +=
            n*Us_[nodei][celli]
           *breakupKernel_->Kb(d, celli)
           *breakupKernel_->massNodeSource  //Birth and death
            (
                abscissa,
                momentOrder
            );
    }

    return cmptMultiply(bSource, validDirections_);
}

void Foam::polydispersePhaseModel::solveSourceOde()
{
    if (!coalescence_ && !breakup_)
    {
        return;
    }

    coalescenceKernel_->preUpdate();
    breakupKernel_->preUpdate();

    volScalarMomentFieldSet& moments = quadrature_.moments();
    label nMoments = quadrature_.nMoments();
    PtrList<volVectorField>& Ups = quadrature_.velocityMoments();
    label nVelocityMoments = Ups.size();

    scalar globalDt = moments[0].mesh().time().deltaT().value();

    if (!solveOde_)
    {
        forAll(moments[0], celli)
        {
            forAll(moments, mi)
            {
                moments[mi][celli] +=
                    globalDt
                   *(
                        coalescenceSource(celli, mi)
                        + breakupSource(celli, mi)
                    );
            }

            forAll(Ups, mi)
            {
                Ups[mi][celli] +=
                    globalDt
                   *(
                        coalescenceSourceU(celli, mi)
                      + breakupSourceU(celli, mi)
                    );
            }
        }
        return;
    }


    Info << "Solving source terms in realizable ODE solver." << endl;

    forAll(moments[0], celli)
    {
        if ((*this)[celli] < 0.01)
        {
            continue;
        }

        scalarField oldMoments(nMoments, Zero);
        vectorField oldUps(nVelocityMoments, Zero);
        forAll(oldMoments, mi)
        {
            oldMoments[mi] = moments[mi][celli];
        }
        forAll(oldUps, mi)
        {
            oldUps[mi] = Ups[mi][celli];
        }

        //- Local time
        scalar localT = 0.0;

        // Initialize the local step
        scalar localDt = localDt_[celli];

        // Initialize RK parameters
        scalarField k1(nMoments_, Zero);
        scalarField k2(nMoments_, Zero);
        scalarField k3(nMoments_, Zero);

        vectorField k1U(nVelocityMoments, Zero);
        vectorField k2U(nVelocityMoments, Zero);
        vectorField k3U(nVelocityMoments, Zero);

        // Flag to indicate if the time step is complete
        bool timeComplete = false;

        // Check realizability of intermediate moment sets
        bool realizableUpdate1 = true;
        bool realizableUpdate2 = true;
        bool realizableUpdate3 = true;

        scalarList momentsSecondStep(nMoments_, 0.0);
        scalar error = 0.0;

        while (!timeComplete)
        {
            bool nullSource = true;
            do
            {
                nullSource = true;

                // First intermediate update
                forAll(oldMoments, mi)
                {
                    k1[mi] =
                        localDt
                       *(
                            coalescenceSource(celli, mi)
                          + breakupSource(celli, mi)
                        );
                    moments[mi][celli] = oldMoments[mi] + k1[mi];
                    nullSource =
                    (
                        mag(k1[mi]) < pow(SMALL, mi)
                     && nullSource
                    );
                }

                if (nullSource)
                {
                    break;
                }

                forAll(oldUps, mi)
                {
                    k1U[mi] =
                        localDt
                       *(
                            coalescenceSourceU(celli, mi)
                          + breakupSourceU(celli, mi)
                        );
                    Ups[mi][celli] = oldUps[mi] + k1U[mi];
                }
                realizableUpdate1 =
                        quadrature_.updateAllLocalQuadrature(celli, false);

                // Second moment update
                forAll(oldMoments, mi)
                {
                    k2[mi] =
                        localDt
                       *(
                            coalescenceSource(celli, mi)
                          + breakupSource(celli, mi)
                        );
                    moments[mi][celli] = oldMoments[mi] + (k1[mi] + k2[mi])/4.0;
                    momentsSecondStep[mi] = moments[mi][celli];
                }
                forAll(oldUps, mi)
                {
                    k2U[mi] =
                        localDt
                       *(
                            coalescenceSourceU(celli, mi)
                          + breakupSourceU(celli, mi)
                        );
                    Ups[mi][celli] = oldUps[mi] + (k1U[mi] + k2U[mi])/4.0;
                }
                realizableUpdate2 =
                        quadrature_.updateAllLocalQuadrature(celli, false);

                // Third moment update
                forAll(oldMoments, mi)
                {
                    k3[mi] =
                        localDt
                       *(
                            coalescenceSource(celli, mi)
                          + breakupSource(celli, mi)
                        );
                    moments[mi][celli] =
                        oldMoments[mi] + (k1[mi] + k2[mi] + 4.0*k3[mi])/6.0;
                }
                forAll(oldUps, mi)
                {
                    k3U[mi] =
                        localDt
                       *(
                            coalescenceSourceU(celli, mi)
                          + breakupSourceU(celli, mi)
                        );
                    Ups[mi][celli] =
                        oldUps[mi] + (k1U[mi] + k2U[mi] + 4.0*k3U[mi])/6.0;
                }
                realizableUpdate3 =
                        quadrature_.updateAllLocalQuadrature(celli, false);

                if
                (
                    !realizableUpdate1
                 || !realizableUpdate2
                 || !realizableUpdate3
                )
                {
                    Info << "Not realizable" << endl;

                    forAll(oldMoments, mi)
                    {
                        moments[mi][celli] = oldMoments[mi];
                    }
                    forAll(oldUps, mi)
                    {
                        Ups[mi][celli] = oldUps[mi];
                    }

                    // Updating local quadrature with old moments
                    quadrature_.updateAllLocalQuadrature(celli);

                    localDt /= 2.0;

                    if (localDt < minLocalDt_)
                    {
                        FatalErrorInFunction
                            << "Reached minimum local step in realizable ODE"
                            << nl
                            << "    solver. Cannot ensure realizability." << nl
                            << abort(FatalError);
                    }
                }
            }
            while
            (
                !realizableUpdate1
             || !realizableUpdate2
             || !realizableUpdate3
            );

            if (nullSource)
            {
                timeComplete = true;
                localT = 0.0;
                break;
            }

            for (label mi = 0; mi < nMoments; mi++)
            {
                scalar scalei =
                   Foam::max
                    (
                        mag(momentsSecondStep[mi]), mag(oldMoments[mi])
                    )*RTol_;

                if (scalei > 0)
                {
                    error +=
                        sqr
                        (
                            (momentsSecondStep[mi] - moments[mi][celli])/scalei
                        );
                }
            }
            error = sqrt(error/nMoments);
            if (error < SMALL)
            {
                timeComplete = true;
                localT = 0.0;
                break;
            }
            else if (error < 1)
            {
                localT += localDt;
                localDt_[celli] = localDt;

                localDt *= Foam::min(facMax_, Foam::max(facMin_, fac_/pow(error, 1.0/3.0)));

                scalar maxLocalDt = Foam::max(globalDt - localT, 0.0);
                localDt = Foam::min(maxLocalDt, localDt);

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }
                forAll(oldUps, mi)
                {
                    oldUps[mi] = Ups[mi][celli];
                }

                if (localDt == 0.0)
                {
                    timeComplete = true;
                    localT = 0.0;
                    break;
                }
            }
            else
            {
                localDt *= Foam::min(1.0, Foam::max(facMin_, fac_/pow(error, 1.0/3.0)));

                forAll(oldMoments, mi)
                {
                    moments[mi][celli] = oldMoments[mi];
                }
                forAll(oldUps, mi)
                {
                    Ups[mi][celli] = oldUps[mi];
                }
                quadrature_.updateAllLocalQuadrature(celli);

                if (localDt < minLocalDt_)
                {
                    Info<<" cell "<<celli
                        <<" error: " <<error
                        <<" dt: "<<localDt<<endl;
                    FatalErrorInFunction
                        << "Reached minimum local step in realizable ODE"
                        << nl
                        << "    solver. Cannot ensure realizability." << nl
                        << abort(FatalError);
                }
            }
        }
    }

    quadrature_.updateAllMoments();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polydispersePhaseModel::polydispersePhaseModel
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
    solveOde_(pbeDict_.lookupOrDefault("ode", false)),
    coalescence_(pbeDict_.lookup("coalescence")),
    breakup_(pbeDict_.lookup("breakup")),
    quadrature_(phaseName, fluid.mesh(), "RPlus"),
    nNodes_(quadrature_.nodes().size()),
    nMoments_(quadrature_.nMoments()),
    alphas_(nNodes_),
    Us_(quadrature_.velocities()),
    Vs_(nNodes_),
    ds_(nNodes_),
    maxD_("maxD", dimLength, phaseDict_),
    minD_("minD", dimLength, phaseDict_),
    minLocalDt_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("minLocalDt"))),
    localDt_(this->size(), fluid.mesh().time().deltaT().value()/10.0),
    ATol_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("RTol"))),
    facMax_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("facMax"))),
    facMin_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("facMin"))),
    fac_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("fac"))),
    validDirections_
    (
        (vector(fluid_.mesh().solutionD()) + vector(1.0, 1.0, 1.0))/2.0
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
    correct();

    const dictionary& pimpleDict =
        fluid_.mesh().solutionDict().subDict("PIMPLE");
    label nCorrectors = pimpleDict.lookupOrDefault<label>("nFluxCorrectors", 0);
    if (nCorrectors > 0)
    {
        word patchName
        (
            pimpleDict.lookupOrDefault
            (
                "corrPatch",
                U_.boundaryField()[0].patch().name()
            )
        );
        wordList boundaries(U_.boundaryField().size(), "zeroGradient");
        forAll(boundaries, patchi)
        {
            if (U_.boundaryField()[patchi].patch().name() == patchName)
            {
                boundaries[patchi] = "fixedFace";
            }
        }

        corr_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "corr",
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid_.mesh(),
                dimensionedScalar("0", sqr(dimLength)/dimTime, 0.0),
                boundaries
            )
        );
    }
}


void Foam::polydispersePhaseModel::setModels()
{
    coalescenceKernel_.set
    (
        new populationBalanceSubModels::aggregationKernels::coalescence
        (
            pbeDict_.subDict("coalescenceKernel"),
            fluid_.mesh()
        )
    );
    breakupKernel_ =
    (
        Foam::populationBalanceSubModels::breakupKernel::New
        (
            pbeDict_.subDict("breakupKernel"),
            fluid_.mesh()
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polydispersePhaseModel::~polydispersePhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polydispersePhaseModel::correct()
{

    if (nNodes_ == 1)
    {
        alphas_[0] = *this;
        ds_[0] =
            Foam::min
            (
                Foam::max
                (
                    Foam::pow
                    (
                        quadrature_.nodes()[0].primaryAbscissae()[0]*6.0
                       /(rho()*Foam::constant::mathematical::pi)
                      + dimensionedScalar("smallVolume", dimVolume, SMALL),
                        1.0/3.0
                    ),
                    minD_
                ),
                maxD_
            );
        d_ = ds_[0];
        return;
    }
    else
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
            alphas_[nodei] =
                node.primaryWeight()*node.primaryAbscissae()[0]/rho()*scale;
            alphas_[nodei].max(0);
            alphas_[nodei].min(1);

            //  Calculate bubble diameter based on bubble mass (abscissa)
            ds_[nodei] =
                Foam::min
                (
                    Foam::max
                    (
                        Foam::pow
                        (
                            node.primaryAbscissae()[0]*6.0
                           /(rho()*Foam::constant::mathematical::pi)
                          + dimensionedScalar("smallVolume", dimVolume, SMALL),
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
}


void Foam::polydispersePhaseModel::relativeTransport()
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
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phiv, zeroPhi)
              + nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
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
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phiv, zeroPhi)
              + quadrature_.velocitiesOwn()[nodei]
               *nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
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
    correct();
}


void Foam::polydispersePhaseModel::averageTransport
(
    const PtrList<fvVectorMatrix>& AEqns
)
{
    // Correct mean flux
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_.nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_.nodesNei();
    dimensionedScalar zeroPhi("zero", phiPtr_().dimensions(), 0.0);
    surfaceScalarField phi(phiPtr_());

    const dictionary& pimpleDict =
        fluid_.mesh().solutionDict().subDict("PIMPLE");
    label nCorrectors = pimpleDict.lookupOrDefault<label>("nFluxCorrectors", 0);
    if (corr_.valid())
    {
        volScalarField& corr = corr_.ref();

        for (label i = 0; i < nCorrectors; i++)
        {
            //- Update interpolated nodes since upwinding direction
            //  may have changed
            quadrature_.interpolateNodes();
            volScalarField meanM1Flux
            (
                IOobject
                (
                    "meanM1Flux",
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, Zero)
            );

            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
                meanM1Flux +=
                    fvc::surfaceIntegrate
                    (
                        nodesNei[nodei].primaryWeight()
                       *nodesNei[nodei].primaryAbscissae()[0]
                       *Foam::min(phi, zeroPhi)
                      + nodesOwn[nodei].primaryWeight()
                       *nodesOwn[nodei].primaryAbscissae()[0]
                       *Foam::max(phi, zeroPhi)
                    );
            }

            dimensionedScalar minCoeff
            (
                dimensionedScalar::lookupOrDefault
                (
                    "minCoeff",
                    pimpleDict,
                    dimDensity,
                    residualAlpha_.value()
                )
            );

            fvScalarMatrix corrEqn
            (
                ((*this)*rho() - quadrature_.moments()[1])
               /fluid_.mesh().time().deltaT()
              + meanM1Flux
              + fvm::laplacian
                (
                    Foam::max(quadrature_.moments()[1], minCoeff),
                    corr,
                    "laplacian(" + quadrature_.moments()[1].name() + ",corr)"
                )
            );
            corrEqn.relax();
            corrEqn.solve();
            phi += fvc::snGrad(corr)*fluid_.mesh().magSf();
        }
    }
    quadrature_.interpolateNodes();

    // Mean moment advection
    Info<< "Transporting moments with average velocity" << endl;
    forAll(quadrature_.moments(), mEqni)
    {
        volScalarField& m = quadrature_.moments()[mEqni];

        volScalarField meanDivUbMp
        (
            IOobject
            (
                "meanDivUbMp",
                fluid_.mesh().time().timeName(),
                fluid_.mesh()
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
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phi, zeroPhi)
              + nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
                    mEqni
                )*Foam::max(phi, zeroPhi)
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

        quadrature_.updateQuadrature();
        Us_[0] = U_;
        return;
    }

    forAll(quadrature_.velocityMoments(), mEqni)
    {
        volVectorField& Up = quadrature_.velocityMoments()[mEqni];

        volVectorField meanDivUbUp
        (
            IOobject
            (
                "meanDivUbUp",
                fluid_.mesh().time().timeName(),
                fluid_.mesh()
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
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phi, zeroPhi)
              + quadrature_.velocitiesOwn()[nodei]
               *nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
                    mEqni
                )*Foam::max(phi, zeroPhi)
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
            Foam::max
            (
                (0.5 + 0.5*tanh(((*this) - 0.63)/0.01))*GREAT,
                residualAlpha_
            )
        );
        tauC.dimensions().reset(dimDensity/dimTime);

        volScalarField alphaRhoi(alphas_[nodei]*rho());

        // Solve for velocities using acceleration terms
        fvVectorMatrix UsEqn
        (
            alphaRhoi*fvm::ddt(Us_[nodei])
          - alphaRhoi*fvc::ddt(Us_[nodei])
          + fvm::Sp(tauC, Us_[nodei])
         ==
            AEqns[nodei]
          + tauC*U_
        );

        UsEqn.relax();
        UsEqn.solve();
    }
    quadrature_.updateAllMoments();

    // Update moments with breakup and coalescence sources
    solveSourceOde();

    //- Update mean velocity
    this->updateVelocity();

    // Update deviation velocity
    forAll(Vs_, nodei)
    {
        Vs_[nodei] = Us_[nodei] - U_;
    }
}


bool Foam::polydispersePhaseModel::read(const bool readOK)
{
    bool read = false;
    if (readOK)
    {
        maxD_.readIfPresent(phaseDict_);
        minD_.readIfPresent(phaseDict_);
        read = true;
    }

    if (pbeDict_.modified())
    {
        const dictionary& odeDict(pbeDict_.subDict("odeCoeffs"));
        pbeDict_.lookup("coalescence") >> coalescence_;
        pbeDict_.lookup("breakup") >> breakup_;
        odeDict.lookup("minLocalDt") >> minLocalDt_;
        odeDict.lookup("ATol") >> ATol_;
        odeDict.lookup("RTol") >> RTol_;
        odeDict.lookup("facMax") >> facMax_;
        odeDict.lookup("facMin") >> facMin_;
        odeDict.lookup("fac") >> fac_;
        read = true;
    }

    return read || readOK;
}
// ************************************************************************* //
