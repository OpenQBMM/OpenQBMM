/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

#include "AGmomentTransportModel.H"
#include "fvc.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"
#include "emptyFvPatch.H"
#include "coupledFvPatch.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AGmomentTransportModel::AGmomentTransportModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const phaseModel& phase,
    volScalarField& Theta,
    volSymmTensorField& Sigma
)
:
    mesh_(mesh),
    phase_(phase),
    alphap_
    (
        mesh.lookupObjectRef<volScalarField>
        (
            IOobject::groupName("alpha", phase.name())
        )
    ),
    Up_
    (
        mesh.lookupObjectRef<volVectorField>
        (
            IOobject::groupName("U", phase.name())
        )
    ),
    Theta_(Theta),
    Sigma_(Sigma),
    Pp_
    (
        IOobject
        (
            IOobject::groupName("Pp", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Theta_*symmTensor::I - Sigma_,
        Sigma_.boundaryField().types()
    ),
    hq_(3, dict.lookupOrDefault<label>("nHerNodePerDim", 4)),
    weights_(hq_.hermiteWeights().size()),
    abscissae_(hq_.hermiteAbscissae().size()),
    ew_(dict.lookupOrDefault<scalar>("wallRestitutionCoefficient", 1)),
    phiw_(dict.lookupOrDefault<scalar>("wallSpecularityCoefficient", 0)),
    own_
    (
        IOobject
        (
            "own",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("own", dimless, 1.0)
    ),
    nei_
    (
        IOobject
        (
            "nei",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("nei", dimless, -1.0)
    ),
    F0_
    (
        IOobject
        (
            "F0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless*dimVol/dimTime, scalar(0))
    ),
    F1_
    (
        IOobject
        (
            "F1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity*dimVol/dimTime, vector::zero)
    ),
    F2_
    (
        IOobject
        (
            "dM2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor
        (
            "zero",
            pow(dimVelocity,2)*dimVol/dimTime,
            symmTensor::zero
        )
    ),
    pDxyz_
    (
        IOobject
        (
            "dxyz",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        vector::zero
    ),
    ddtAlphaDilute_
    (
        IOobject
        (
            "ddtAlphaDilute",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", inv(dimTime), 0.0)
    )
{
    forAll(weights_, nodei)
    {
        weights_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "weight" + Foam::name(nodei),
                        phase.name()
                    ),
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        abscissae_.set
        (
            nodei,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "abscisa" + Foam::name(nodei),
                        phase.name()
                    ),
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector("zero", dimVelocity, Zero)
            )
        );
    }
    const faceList & ff = mesh.faces();
    const pointField & pp = mesh.points();

    forAll ( mesh.C(), celli)
    {
        const cell & cc = mesh_.cells()[celli];
        labelList pLabels(cc.labels(ff));
        pointField pLocal(pLabels.size(), vector::zero);

        forAll(pLabels, pointi)
        {
            pLocal[pointi] = pp[pLabels[pointi]];
        }

        pDxyz_[celli][0] =
            Foam::max(pLocal & vector(1,0,0))
          - Foam::min(pLocal & vector(1,0,0));

        pDxyz_[celli][1] =
            Foam::max(pLocal & vector(0,1,0))
          - Foam::min(pLocal & vector(0,1,0));

        pDxyz_[celli][2] =
            Foam::max(pLocal & vector(0,0,1))
          - Foam::min(pLocal & vector(0,0,1));
    }

    forAll(mesh_.boundary(), patchi)
    {
        pDxyz_.boundaryFieldRef()[patchi] == vector(1,1,1);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::AGmomentTransportModel::~AGmomentTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar  Foam::AGmomentTransportModel::maxUxDx() const
{
    scalar maxUxDx = 0.0;
    forAll(weights_, nodei)
    {
        maxUxDx = max(maxUxDx, max(mag(abscissae_[nodei])).value());
    }

    return maxUxDx;
}


void Foam::AGmomentTransportModel::solve
(
    const surfaceScalarField& h2f
)
{
    tmp<surfaceScalarField> th1f(1.0 - h2f);
    const surfaceScalarField& h1f = th1f();

    Pp_ = Theta_*symmTensor::I - Sigma_;
    Pp_.correctBoundaryConditions();

    calcMomentFluxes(h1f);

    volScalarField m0
    (
        IOobject::groupName("moment.0", phase_.name()),
        alphap_
    );

    volVectorField m1
    (
        IOobject::groupName("moment.1", phase_.name()),
        alphap_*Up_
    );

    volSymmTensorField m2
    (
        IOobject::groupName("moment.2", phase_.name()),
        alphap_*(Pp_ + sqr(Up_))
    );

    volScalarField m0Old(m0);
    volVectorField m1Old(m1);
    volSymmTensorField m2Old(m2);

    const dimensionedScalar& deltaT = mesh_.time().deltaT();

    // Predictor step
    m0 = m0Old - 0.5*fvc::surfaceIntegrate(F0_)*deltaT;
    m1 = m1Old - 0.5*fvc::surfaceIntegrate(F1_)*deltaT;
    m2 = m2Old - 0.5*fvc::surfaceIntegrate(F2_)*deltaT;

    alphap_ = m0;
    alphap_.correctBoundaryConditions();

    Up_ = m1/max(m0, phase_.residualAlpha());
    Up_.correctBoundaryConditions();

    Pp_ = m2/max(m0, phase_.residualAlpha()) - sqr(Up_);
    Pp_.correctBoundaryConditions();

    calcMomentFluxes(h1f);

    // Correction
    m0 = m0Old - fvc::surfaceIntegrate(F0_)*deltaT;
    m1 = m1Old - fvc::surfaceIntegrate(F1_)*deltaT;
    m2 = m2Old - fvc::surfaceIntegrate(F2_)*deltaT;

    // Set volume fraction updated form dilute transport
    alphap_ = m0;
    alphap_.correctBoundaryConditions();
    ddtAlphaDilute_ = fvc::ddt(alphap_);
    alphap_.storeOldTime();

    // Set velocity from dilute transport
    Up_ = m1/max(m0, phase_.residualAlpha());
    Up_.correctBoundaryConditions();
    Up_.storeOldTime();

    // Update particle pressure tensor
    Pp_ = m2/max(m0, phase_.residualAlpha()) - sqr(Up_);

    forAll(Pp_,i)
    {
        if (Pp_[i].xx() < SMALL)
        {
            Pp_[i].xx() = SMALL;
        }

        if (Pp_[i].yy() < SMALL)
        {
            Pp_[i].yy() = SMALL;
        }

        if (Pp_[i].zz() < SMALL)
        {
            Pp_[i].zz() = SMALL;
        }
    }

    Pp_.correctBoundaryConditions();

    // Update granular temperature based on granular temperature
    Theta_ = 1.0/3.0*tr(Pp_);
    Theta_.max(0);
    Theta_.min(100);
    Theta_.correctBoundaryConditions();
    Theta_.storeOldTime();

    // Update granular stress tensor
    Sigma_ = Theta_*symmTensor::I - Pp_;
    Sigma_.correctBoundaryConditions();
    Sigma_.storeOldTime();
}


void Foam::AGmomentTransportModel::calcMomentFluxes
(
    const surfaceScalarField& h1f
)
{
    PtrList<surfaceScalarField> weightsOwn(weights_.size());
    PtrList<surfaceScalarField> weightsNei(weights_.size());
    PtrList<surfaceVectorField> abscissaeOwn(abscissae_.size());
    PtrList<surfaceVectorField> abscissaeNei(abscissae_.size());

    forAll(alphap_, celli)
    {
        hq_.calcHermiteQuadrature(Up_[celli], Pp_[celli]);
        forAll(weights_, nodei)
        {
            weights_[nodei][celli] = hq_.hermiteWeights()[nodei]*alphap_[celli];
            abscissae_[nodei][celli] = hq_.hermiteAbscissae()[nodei];
        }
    }
    forAll(alphap_.boundaryField(), patchi)
    {
        forAll(alphap_.boundaryField()[patchi], facei)
        {
            hq_.calcHermiteQuadrature
            (
                Up_.boundaryField()[patchi][facei],
                Pp_.boundaryField()[patchi][facei]
            );
            forAll(weights_, nodei)
            {
                weights_[nodei].boundaryFieldRef()[patchi][facei] =
                    hq_.hermiteWeights()[nodei]*alphap_.boundaryField()[patchi][facei];
                abscissae_[nodei].boundaryFieldRef()[patchi][facei] =
                    hq_.hermiteAbscissae()[nodei];
            }
        }
    }

    forAll(weights_, nodei)
    {
        weights_[nodei].correctBoundaryConditions();
        abscissae_[nodei].correctBoundaryConditions();
    }

    forAll(weightsOwn, nodei)
    {
        weightsOwn.set
        (
            nodei,
            fvc::interpolate(weights_[nodei], own_, "reconstruct(weight)")
        );
        weightsNei.set
        (
            nodei,
            fvc::interpolate(weights_[nodei], nei_, "reconstruct(weight)")
        );

        abscissaeOwn.set
        (
            nodei,
            fvc::interpolate(abscissae_[nodei], own_, "reconstruct(U)")
        );
        abscissaeNei.set
        (
            nodei,
            fvc::interpolate(abscissae_[nodei], nei_, "reconstruct(U)")
        );
    }

    // Update face values for wall conditions
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& currPatch = mesh_.boundary()[patchi];
        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField& bfSf(mesh_.Sf().boundaryField()[patchi]);
            vectorField bfNorm(bfSf/mag(bfSf));

            forAll(weights_, nodei)
            {

                const volScalarField& weight = weights_[nodei];
                surfaceScalarField& weightOwn = weightsOwn[nodei];
                surfaceScalarField& weightNei = weightsNei[nodei];
                const volVectorField& U = abscissae_[nodei];
                surfaceVectorField& UOwn = abscissaeOwn[nodei];;
                surfaceVectorField& UNei = abscissaeNei[nodei];

                scalarField& bfwOwn = weightOwn.boundaryFieldRef()[patchi];
                scalarField& bfwNei = weightNei.boundaryFieldRef()[patchi];
                vectorField& bfUOwn = UOwn.boundaryFieldRef()[patchi];
                vectorField& bfUNei = UNei.boundaryFieldRef()[patchi];

                forAll(currPatch, facei)
                {
                    label faceCelli = currPatch.faceCells()[facei];

                    bfwOwn[facei] = weight[faceCelli];
                    bfUOwn[facei] = U[faceCelli];

                    bfwNei[facei] = weight[faceCelli];
                    bfUNei[facei] =
                        U[faceCelli]
                      - (1.0 + ew_)*(U[faceCelli] & bfNorm[facei])
                       *bfNorm[facei];
                }
            }
        }
    }

    dimensionedScalar zeroPhi("zeroPhi", dimVelocity*dimArea, 0.0);
    F0_ == dimensionedScalar("zero", F0_.dimensions(), 0.0);
    F1_ == dimensionedVector("zero", F1_.dimensions(), Zero);
    F2_ == dimensionedSymmTensor("zero", F2_.dimensions(), Zero);

    forAll(weights_, nodei)
    {
        surfaceScalarField phiOwn(abscissaeOwn[nodei] & mesh_.Sf());
        surfaceScalarField phiNei(abscissaeNei[nodei] & mesh_.Sf());

        F0_ +=
            weightsOwn[nodei]*max(phiOwn, zeroPhi)
          + weightsNei[nodei]*min(phiNei, zeroPhi);

        F1_ +=
            weightsOwn[nodei]*abscissaeOwn[nodei]*max(phiOwn, zeroPhi)
          + weightsNei[nodei]*abscissaeNei[nodei]*min(phiNei, zeroPhi);

        F2_ +=
            weightsOwn[nodei]*max(phiOwn, zeroPhi)*sqr(abscissaeOwn[nodei])
          + weightsNei[nodei]*min(phiNei, zeroPhi)*sqr(abscissaeNei[nodei]);
    }

    F0_ *= h1f;
    F1_ *= h1f;
    F2_ *= h1f;

    return;
}

// ************************************************************************* //
