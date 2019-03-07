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

#include "velocityMomentAdvection.H"
#include "IOmanip.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityMomentAdvection, 0);
    defineRunTimeSelectionTable(velocityMomentAdvection, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityMomentAdvection::velocityMomentAdvection
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const word& support
)
:
    name_(quadrature.name()),
    moments_(quadrature.moments()),
    nMoments_(moments_.size()),
    own_
    (
        IOobject
        (
            "own",
            moments_(0).mesh().time().timeName(),
            moments_(0).mesh()
        ),
        moments_(0).mesh(),
        dimensionedScalar("own", dimless, 1.0)
    ),
    nei_
    (
        IOobject
        (
            "nei",
            moments_(0).mesh().time().timeName(),
            moments_(0).mesh()
        ),
        moments_(0).mesh(),
        dimensionedScalar("nei", dimless, -1.0)
    ),
    support_(support),
    momentOrders_(quadrature.momentOrders()),
    nodeIndexes_(quadrature.nodeIndexes()),
    divMoments_(nMoments_),
    ew_(dict.lookupOrDefault("ew", 1.0))
{
    forAll(divMoments_, momenti)
    {
        const labelList& momentOrder = momentOrders_[momenti];
        divMoments_.set
        (
            momenti,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "divMoment"
                       + mappedList<vector>::listToWord(momentOrder),
                        name_
                    ),
                    moments_(0).mesh().time().timeName(),
                    moments_(0).mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                moments_(0).mesh(),
                dimensionedScalar
                (
                    "zero", moments_[momenti].dimensions()/dimTime, 0
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityMomentAdvection::~velocityMomentAdvection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityMomentAdvection::updateWallCollisions
(
    const PtrList<volVelocityNode>& nodes,
    PtrList<surfaceVelocityNode>& nodesOwn,
    PtrList<surfaceVelocityNode>& nodesNei
)
{
    const fvMesh& mesh = own_.mesh();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField& bfSf(mesh.Sf().boundaryField()[patchi]);
            vectorField bfNorm(bfSf/mag(bfSf));

            forAll(nodes, nodei)
            {
                const volVelocityNode& node = nodes[nodei];
                surfaceVelocityNode& nodeNei(nodesNei[nodei]);
                surfaceVelocityNode& nodeOwn(nodesOwn[nodei]);

                const volScalarField& weight = node.primaryWeight();
                surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
                surfaceScalarField& weightNei = nodeNei.primaryWeight();
                const volVectorField& U = node.velocityAbscissae();
                surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
                surfaceVectorField& UNei = nodeNei.velocityAbscissae();

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
                      - (1.0 + this->ew_)*(bfUOwn[facei] & bfNorm[facei])
                       *bfNorm[facei];
                }
            }
        }
    }
}


// ************************************************************************* //
