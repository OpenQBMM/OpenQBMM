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
    PtrList<dimensionSet> abscissaeDimensions(momentOrders_[0].size());
    labelList zeroOrder(momentOrders_[0].size(), 0);

    forAll(abscissaeDimensions, dimi)
    {
        labelList firstOrder(zeroOrder);
        firstOrder[dimi] = 1;

        abscissaeDimensions.set
        (
            dimi,
            new dimensionSet
            (
                moments_(firstOrder).dimensions()/moments_(0).dimensions()
            )
        );
    }

    nodesNei_ = autoPtr<PtrList<surfaceVelocityNode> >
    (
        new PtrList<surfaceVelocityNode>(nodes_.size())
    );

    nodesOwn_ = autoPtr<PtrList<surfaceVelocityNode> >
    (
        new PtrList<surfaceVelocityNode>(nodes_.size())
    );

    PtrList<surfaceVelocityNode>& nodesNei = nodesNei_();
    PtrList<surfaceVelocityNode>& nodesOwn = nodesOwn_();

    // Populating nodes and interpolated nodes
    forAll(nodes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        nodesNei.set
        (
            nodei,
            new surfaceVelocityNode
            (
                "nodeNei" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_(0).mesh(),
                moments_(0).dimensions(),
                abscissaeDimensions,
                false
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceVelocityNode
            (
                "nodeOwn" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_(0).mesh(),
                moments_(0).dimensions(),
                abscissaeDimensions,
                false
            )
        );
    }

    {
        IStringStream weightLimiter("upwind");
        IStringStream scalarAbscissaeLimiter("upwind");
        IStringStream velocityAbscissaeLimiter("upwind");
        weightOwnScheme_ = fvc::scheme<scalar>(own_, weightLimiter);
        scalarAbscissaeOwnScheme_ =
            fvc::scheme<scalar>
            (
                own_,
                scalarAbscissaeLimiter
            );
        velocityAbscissaeOwnScheme_ =
            fvc::scheme<vector>
            (
                own_,
                velocityAbscissaeLimiter
            );
    }

    {
        IStringStream weightLimiter("upwind");
        IStringStream scalarAbscissaeLimiter("upwind");
        IStringStream velocityAbscissaeLimiter("upwind");
        weightNeiScheme_ = fvc::scheme<scalar>(nei_, weightLimiter);
        scalarAbscissaeNeiScheme_ =
            fvc::scheme<scalar>
            (
                nei_,
                scalarAbscissaeLimiter
            );
        velocityAbscissaeNeiScheme_ =
            fvc::scheme<vector>
            (
                nei_,
                velocityAbscissaeLimiter
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityAdvection::firstOrderKinetic::~firstOrderKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityAdvection::firstOrderKinetic::interpolateNodes()
{
    PtrList<surfaceVelocityNode>& nodesNei = nodesNei_();
    PtrList<surfaceVelocityNode>& nodesOwn = nodesOwn_();

    forAll(nodes_, nodei)
    {
        const volVelocityNode& node(nodes_[nodei]);
        surfaceVelocityNode& nodeNei(nodesNei[nodei]);
        surfaceVelocityNode& nodeOwn(nodesOwn[nodei]);

        nodeOwn.primaryWeight() =
            weightOwnScheme_().interpolate(node.primaryWeight());
        nodeOwn.velocityAbscissae() =
            velocityAbscissaeOwnScheme_().interpolate(node.velocityAbscissae());

        nodeNei.primaryWeight() =
            weightNeiScheme_().interpolate(node.primaryWeight());
        nodeNei.velocityAbscissae() =
            velocityAbscissaeNeiScheme_().interpolate(node.velocityAbscissae());

        forAll(node.primaryAbscissae(), cmpt)
        {
            nodeOwn.primaryAbscissae()[cmpt] =
                scalarAbscissaeOwnScheme_().interpolate
                (
                    node.primaryAbscissae()[cmpt]
                );
            nodeNei.primaryAbscissae()[cmpt] =
                scalarAbscissaeNeiScheme_().interpolate
                (
                    node.primaryAbscissae()[cmpt]
                );
        }
    }
}


Foam::scalar
Foam::velocityAdvection::firstOrderKinetic::realizableCo() const
{
    const fvMesh& mesh = this->own_.mesh();
    surfaceVectorField Sf(mesh.Sf());

    scalarField maxCoNum(mesh.nCells(), 1.0);

    forAll(this->nodes_, nodei)
    {

        surfaceScalarField phiOwn
        (
            mag(this->nodesOwn_()[nodei].velocityAbscissae() & mesh.Sf())
        );
        surfaceScalarField phiNei
        (
            mag(this->nodesNei_()[nodei].velocityAbscissae() & mesh.Sf())
        );

        forAll(moments_[0], celli)
        {
            const labelList& cell = mesh.cells()[celli];

            scalar den = 0;
            forAll(cell, facei)
            {
                if (cell[facei] < mesh.nInternalFaces())
                {
                    den +=
                        max
                        (
                            phiOwn[cell[facei]],
                            phiNei[cell[facei]]
                        );
                }

                den = max(den, small);
                maxCoNum[celli] =
                    min
                    (
                        maxCoNum[celli],
                        mesh.V()[celli]
                       /(den*mesh.time().deltaTValue())
                    );
            }
        }
    }
    return gMin(maxCoNum);
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
                        mag(fvc::flux(nodes_[nodei].velocityAbscissae()))
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
    updateWallCollisions(nodes_, nodesOwn_(), nodesNei_());

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

    const labelList& scalarIndexes = nodes_[0].scalarIndexes();
    const labelList& velocityIndexes = nodes_[0].velocityIndexes();

    forAll(nodes_, nodei)
    {
        const surfaceVelocityNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVelocityNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();

        const PtrList<surfaceScalarField>& scalarAbscissaeOwn =
            nodeOwn.primaryAbscissae();
        const PtrList<surfaceScalarField>& scalarAbscissaeNei =
            nodeNei.primaryAbscissae();

        const surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        const surfaceVectorField& UNei = nodeNei.velocityAbscissae();

        surfaceScalarField phiOwn(UOwn & mesh.Sf());
        surfaceScalarField phiNei(UNei & mesh.Sf());

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

            forAll(scalarIndexes, cmpti)
            {
                const label cmpt = scalarIndexes[cmpti];
                const label cmptMomentOrder = momentOrder[cmpt];

                if (cmptMomentOrder > 0)
                {
                    const surfaceScalarField& abscissaOwnCmpt =
                    scalarAbscissaeOwn[cmpti];
                    const surfaceScalarField& abscissaNeiCmpt =
                        scalarAbscissaeNei[cmpti];

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
            }

            forAll(velocityIndexes, cmpti)
            {
                const label cmpt = velocityIndexes[cmpti];
                const label cmptMomentOrder = momentOrder[cmpt];

                if (cmptMomentOrder > 0)
                {
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
        updateWallCollisions(nodes_, nodesOwn_(), nodesNei_());
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

    const labelList& scalarIndexes = nodes_[0].scalarIndexes();
    const labelList& velocityIndexes = nodes_[0].velocityIndexes();

    forAll(nodes_, nodei)
    {
        const surfaceVelocityNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVelocityNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();

        const PtrList<surfaceScalarField>& scalarAbscissaeOwn =
            nodeOwn.primaryAbscissae();
        const PtrList<surfaceScalarField>& scalarAbscissaeNei =
            nodeNei.primaryAbscissae();

        const surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        const surfaceVectorField& UNei = nodeNei.velocityAbscissae();

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

            forAll(scalarIndexes, cmpti)
            {
                const label cmpt = scalarIndexes[cmpti];
                const label cmptMomentOrder = momentOrder[cmpt];

                if (cmptMomentOrder > 0)
                {
                    const surfaceScalarField& abscissaOwnCmpt =
                    scalarAbscissaeOwn[cmpti];
                    const surfaceScalarField& abscissaNeiCmpt =
                        scalarAbscissaeNei[cmpti];

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
            }

            forAll(velocityIndexes, cmpti)
            {
                const label cmpt = velocityIndexes[cmpti];
                const label cmptMomentOrder = momentOrder[cmpt];

                if (cmptMomentOrder > 0)
                {
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
        updateWallCollisions(nodes_, nodesOwn_(), nodesNei_());
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

    const labelList& scalarIndexes = nodes_[0].scalarIndexes();
    const labelList& velocityIndexes = nodes_[0].velocityIndexes();

    forAll(nodes_, nodei)
    {
        const surfaceVelocityNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVelocityNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();

        const PtrList<surfaceScalarField>& scalarAbscissaeOwn =
            nodeOwn.primaryAbscissae();
        const PtrList<surfaceScalarField>& scalarAbscissaeNei =
            nodeNei.primaryAbscissae();

        const surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        const surfaceVectorField& UNei = nodeNei.velocityAbscissae();

        surfaceVectorField VOwn
        (
            velocityAbscissaeOwnScheme_().interpolate(Us[nodei])
        );
        surfaceVectorField VNei
        (
            velocityAbscissaeNeiScheme_().interpolate(Us[nodei])
        );
        surfaceScalarField phiOwn(VOwn & mesh.Sf());
        surfaceScalarField phiNei(VNei & mesh.Sf());

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

            forAll(scalarIndexes, cmpti)
            {
                const label cmpt = scalarIndexes[cmpti];
                const label cmptMomentOrder = momentOrder[cmpt];

                if (cmptMomentOrder > 0)
                {
                    const surfaceScalarField& abscissaOwnCmpt =
                    scalarAbscissaeOwn[cmpti];
                    const surfaceScalarField& abscissaNeiCmpt =
                        scalarAbscissaeNei[cmpti];

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
            }

            forAll(velocityIndexes, cmpti)
            {
                const label cmpt = velocityIndexes[cmpti];
                const label cmptMomentOrder = momentOrder[cmpt];

                if (cmptMomentOrder > 0)
                {
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
