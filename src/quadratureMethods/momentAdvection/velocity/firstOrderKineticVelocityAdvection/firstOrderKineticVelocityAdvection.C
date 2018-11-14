/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
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

#include "firstOrderKineticVelocityAdvection.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{
    defineTypeNameAndDebug(firstOrderKinetic, 0);

    addToRunTimeSelectionTable
    (
        velocityMomentAdvection,
        firstOrderKinetic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityAdvection::firstOrderKinetic::firstOrderKinetic
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const word& support
)
:
    velocityMomentAdvection(dict, quadrature, support),
    nodes_(quadrature.nodes()),
    nodesNei_(),
    nodesOwn_()
{
    nodesNei_ = autoPtr<PtrList<surfaceVectorNode> >
    (
        new PtrList<surfaceVectorNode>(nodes_.size())
    );

    nodesOwn_ = autoPtr<PtrList<surfaceVectorNode> >
    (
        new PtrList<surfaceVectorNode>(nodes_.size())
    );

    PtrList<surfaceVectorNode>& nodesNei = nodesNei_();
    PtrList<surfaceVectorNode>& nodesOwn = nodesOwn_();

    // Populating nodes and interpolated nodes
    forAll(nodes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        nodesNei.set
        (
            nodei,
            new surfaceVectorNode
            (
                "nodeNei" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                false
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceVectorNode
            (
                "nodeOwn" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                false
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityAdvection::firstOrderKinetic::~firstOrderKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityAdvection::firstOrderKinetic::interpolateNodes()
{
    PtrList<surfaceVectorNode>& nodesNei = nodesNei_();
    PtrList<surfaceVectorNode>& nodesOwn = nodesOwn_();

    forAll(nodes_, nodei)
    {
        const volVectorNode& node(nodes_[nodei]);
        surfaceVectorNode& nodeNei(nodesNei[nodei]);
        surfaceVectorNode& nodeOwn(nodesOwn[nodei]);

        nodeOwn.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), own_, "reconstruct(weight)");

        nodeOwn.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                own_,
                "reconstruct(U)"
            );

        nodeNei.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), nei_, "reconstruct(weight)");

        nodeNei.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                nei_,
                "reconstruct(U)"
            );
    }
}

void Foam::velocityAdvection::firstOrderKinetic::updateWallCollisions()
{
    const fvMesh& mesh = own_.mesh();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField& bfSf(mesh.Sf().boundaryField()[patchi]);
            vectorField bfNorm(bfSf/mag(bfSf));

            forAll(nodes_, nodei)
            {
                const volVectorNode& node = nodes_[nodei];
                surfaceVectorNode& nodeNei(nodesNei_()[nodei]);
                surfaceVectorNode& nodeOwn(nodesOwn_()[nodei]);

                const volScalarField& weight = node.primaryWeight();
                surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
                surfaceScalarField& weightNei = nodeNei.primaryWeight();
                const volVectorField& U = node.primaryAbscissa();
                surfaceVectorField& UOwn = nodeOwn.primaryAbscissa();
                surfaceVectorField& UNei = nodeNei.primaryAbscissa();

                scalarField& bfwOwn = weightOwn.boundaryFieldRef()[patchi];
                scalarField& bfwNei = weightNei.boundaryFieldRef()[patchi];
                vectorField& bfUOwn = UOwn.boundaryFieldRef()[patchi];
                vectorField& bfUNei = UNei.boundaryFieldRef()[patchi];

                forAll(currPatch, facei)
                {
                    label faceCelli = currPatch.faceCells()[facei];

                    bfwOwn[facei] = weight[faceCelli];
                    bfUOwn[facei] = U[faceCelli];

                    bfwNei[facei] = bfwOwn[facei];
                    bfUNei[facei] =
                        bfUOwn[facei]
                      - 2.0*(bfUOwn[facei] & bfNorm[facei])
                       *bfNorm[facei];
                }
            }
        }
    }
}


Foam::scalar
Foam::velocityAdvection::firstOrderKinetic::realizableCo() const
{
    return 1.0;
}

Foam::scalar Foam::velocityAdvection::firstOrderKinetic::CoNum() const
{
    scalar CoNum = 0.0;
    const fvMesh& mesh = own_.mesh();
    forAll(nodes_, nodei)
    {
        CoNum =
            max
            (
                CoNum,
                0.5*gMax
                (
                    fvc::surfaceSum
                    (
                        mag(fvc::flux(nodes_[nodei].primaryAbscissa()))
                    )().primitiveField()/mesh.V().field()
                )*mesh.time().deltaTValue()
            );
    }
    return CoNum;
}

void Foam::velocityAdvection::firstOrderKinetic::update()
{
    const fvMesh& mesh = own_.mesh();
    dimensionedScalar zeroPhi("zero", dimVolume/dimTime, 0.0);

    // Interpolate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    updateWallCollisions();

    // Zero moment flux
    forAll(divMoments_, divi)
    {
        divMoments_[divi] =
            dimensionedScalar
            (
                "0",
                moments_[divi].dimensions()/dimTime,
                0.0
            );
    }

    forAll(nodes_, nodei)
    {
        const surfaceVectorNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVectorNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();
        const surfaceVectorField& UOwn = nodeOwn.primaryAbscissa();
        const surfaceVectorField& UNei = nodeNei.primaryAbscissa();

        surfaceScalarField phiOwn(UOwn & mesh.Sf());
        surfaceScalarField phiNei(UNei & mesh.Sf());

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

            forAll(momentOrder, cmpti)
            {
                const label cmptMomentOrder = momentOrder[cmpti];

                tmp<surfaceScalarField> abscissaOwnCmpt =
                   UOwn.component(cmpti);
                tmp<surfaceScalarField> abscissaNeiCmpt =
                    UNei.component(cmpti);

                tmp<surfaceScalarField> mOwnPow =
                    momentCmptOwn
                   *pow
                    (
                        abscissaOwnCmpt,
                        cmptMomentOrder
                    );
                tmp<surfaceScalarField> mNeiPow =
                    momentCmptNei
                   *pow
                    (
                        abscissaNeiCmpt,
                        cmptMomentOrder
                    );
                momentCmptOwn.dimensions().reset(mOwnPow().dimensions());
                momentCmptOwn == mOwnPow;

                momentCmptNei.dimensions().reset(mNeiPow().dimensions());
                momentCmptNei == mNeiPow;
            }

            divMoments_[divi] +=
                fvc::surfaceIntegrate
                (
                    momentCmptOwn*max(phiOwn, zeroPhi)
                  + momentCmptNei*min(phiNei, zeroPhi)
                );
        }
    }
}

void Foam::velocityAdvection::firstOrderKinetic::update
(
    const surfaceScalarField& phi,
    const bool wallCollisions
)
{
    dimensionedScalar zeroPhi("zero", dimVolume/dimTime, 0.0);

    // Interpolate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    if (wallCollisions)
    {
        updateWallCollisions();
    }

    // Zero moment fluxes
    forAll(divMoments_, divi)
    {
        divMoments_[divi] =
            dimensionedScalar
            (
                "0",
                moments_[divi].dimensions()/dimTime,
                0.0
            );
    }

    forAll(nodes_, nodei)
    {
        const surfaceVectorNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVectorNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();
        const surfaceVectorField& abscissaOwn = nodeOwn.primaryAbscissa();
        const surfaceVectorField& abscissaNei = nodeNei.primaryAbscissa();

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

            forAll(momentOrder, cmpti)
            {
                const label cmptMomentOrder = momentOrder[cmpti];

                tmp<surfaceScalarField> abscissaOwnCmpt =
                   abscissaOwn.component(cmpti);
                tmp<surfaceScalarField> abscissaNeiCmpt =
                    abscissaNei.component(cmpti);

                tmp<surfaceScalarField> mOwnPow =
                    momentCmptOwn
                   *pow
                    (
                        abscissaOwnCmpt,
                        cmptMomentOrder
                    );
                tmp<surfaceScalarField> mNeiPow =
                    momentCmptNei
                   *pow
                    (
                        abscissaNeiCmpt,
                        cmptMomentOrder
                    );
                momentCmptOwn.dimensions().reset(mOwnPow().dimensions());
                momentCmptOwn == mOwnPow;

                momentCmptNei.dimensions().reset(mNeiPow().dimensions());
                momentCmptNei == mNeiPow;
            }

            divMoments_[divi] +=
                fvc::surfaceIntegrate
                (
                    momentCmptOwn*max(phi, zeroPhi)
                  + momentCmptNei*min(phi, zeroPhi)
                );
        }
    }
}

void Foam::velocityAdvection::firstOrderKinetic::update
(
    const mappedPtrList<volVectorField>& Us,
    const bool wallCollisions
)
{
    const fvMesh& mesh = own_.mesh();
    dimensionedScalar zeroPhi("zero", dimVolume/dimTime, 0.0);

    // Interplate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    if (wallCollisions)
    {
        updateWallCollisions();
    }

    // Zero moment fluxes
    forAll(divMoments_, divi)
    {
        divMoments_[divi] =
            dimensionedScalar
            (
                "0",
                moments_[divi].dimensions()/dimTime,
                0.0
            );
    }

    forAll(nodes_, nodei)
    {
//         const labelList& nodeIndex = nodeIndexes_[nodei];

        const surfaceVectorNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVectorNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();
        const surfaceVectorField& abscissaOwn = nodeOwn.primaryAbscissa();
        const surfaceVectorField& abscissaNei = nodeNei.primaryAbscissa();

        surfaceVectorField UOwn
        (
            fvc::interpolate(Us[nodei], own_, "reconstruct(U)")
        );
        surfaceVectorField UNei
        (
            fvc::interpolate(Us[nodei], nei_, "reconstruct(U)")
        );

        surfaceScalarField phiOwn(UOwn & mesh.Sf());
        surfaceScalarField phiNei(UNei & mesh.Sf());

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];
            // Calculate size moment flux
            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

            forAll(momentOrder, cmpti)
            {
                const label cmptMomentOrder = momentOrder[cmpti];

                tmp<surfaceScalarField> abscissaOwnCmpt =
                   abscissaOwn.component(cmpti);
                tmp<surfaceScalarField> abscissaNeiCmpt =
                    abscissaNei.component(cmpti);

                tmp<surfaceScalarField> mOwnPow =
                    momentCmptOwn
                   *pow
                    (
                        abscissaOwnCmpt,
                        cmptMomentOrder
                    );
                tmp<surfaceScalarField> mNeiPow =
                    momentCmptNei
                   *pow
                    (
                        abscissaNeiCmpt,
                        cmptMomentOrder
                    );
                momentCmptOwn.dimensions().reset(mOwnPow().dimensions());
                momentCmptOwn == mOwnPow;

                momentCmptNei.dimensions().reset(mNeiPow().dimensions());
                momentCmptNei == mNeiPow;
            }

            divMoments_[divi] +=
                fvc::surfaceIntegrate
                (
                    momentCmptOwn*max(phiOwn, zeroPhi)
                  + momentCmptNei*min(phiNei, zeroPhi)
                );
        }
    }
}

// ************************************************************************* //
