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
    fixedWalls_(moments_[0].boundaryField().size(), false),
    wallTemperatures_(fixedWalls_.size(), 0.0),
    ew_(dict.lookupOrDefault("ew", 1.0))
{
    if (dict.found("fixedTemperatureBoundaries"))
    {
        wordList tmpWalls(dict.lookup("fixedTemperatureBoundaries"));
        scalarList tmpWallTemperatures(dict.lookup("wallTemperatures"));

        const fvMesh& mesh = moments_[0].mesh();
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& currPatch = mesh.boundary()[patchi];
            forAll(tmpWalls, walli)
            {
                if (tmpWalls[walli] == currPatch.name())
                {
                    fixedWalls_[patchi] = true;
                    wallTemperatures_[patchi] = tmpWallTemperatures[walli];
                }
            }
        }
    }


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
    const volScalarField* ThetaPtr;
    word ThetaName(IOobject::groupName("Theta", moments_[0].group()));
    if (mesh.foundObject<volScalarField>(ThetaName))
    {
        ThetaPtr =
        (
            &mesh.lookupObject<volScalarField>(ThetaName)
        );
    }
    const volScalarField& Theta = *ThetaPtr;

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField& bfSf(mesh.Sf().boundaryField()[patchi]);
            vectorField bfNorm(bfSf/mag(bfSf));
            scalarField scale(bfSf.size(), 1.0);

            if (fixedWalls_[patchi] && ThetaPtr)
            {
                scale =
                (
                    sqrt
                    (
                        wallTemperatures_[patchi]
                       /max(Theta.boundaryField()[patchi], 1e-8)
                    )
                );
            }

            scalarField Gin(bfSf.size(), 0.0);
            scalarField Gout(bfSf.size(), 0.0);

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

                bfwOwn = weight.boundaryField()[patchi].patchInternalField();
                bfwNei = bfwOwn;

                bfUOwn = U.boundaryField()[patchi].patchInternalField();
                bfUNei =
                    (
                        bfUOwn
                      - (1.0 + this->ew_)*(bfUOwn & bfNorm)*bfNorm
                    );
                vectorField vn((bfUNei & bfNorm)*bfNorm);
                vectorField vt(bfUNei - vn);
                bfUNei = vn*scale + vt;

                Gin += max(0.0, bfUOwn & bfSf)*bfwOwn;
                Gout -= min(0.0, bfUNei & bfSf)*bfwNei;
            }

            forAll(nodes, nodei)
            {
                surfaceVelocityNode& nodeNei(nodesNei[nodei]);

                scalarField& bfWNei =
                    nodeNei.primaryWeight().boundaryFieldRef()[patchi];

                bfWNei *= Gin/max(Gout, small);
            }
        }
    }
}


// ************************************************************************* //
