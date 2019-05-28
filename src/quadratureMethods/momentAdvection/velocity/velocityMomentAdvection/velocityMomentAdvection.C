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
#include "symmetryFvPatch.H"

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
    labelList velocityIndexes = nodes[0].velocityIndexes();
    label nd = velocityIndexes.size();

    labelList order000(momentOrders_[0].size(), 0);
    labelList order100(order000);
    labelList order010(order000);
    labelList order001(order000);
    labelList order200(order000);
    labelList order020(order000);
    labelList order002(order000);

    order100[velocityIndexes[0]] = 1;
    order200[velocityIndexes[0]] = 2;
    if (nd > 1)
    {
        order010[velocityIndexes[1]] = 1;
        order020[velocityIndexes[1]] = 2;
    }
    if (nd > 2)
    {
        order001[velocityIndexes[2]] = 1;
        order002[velocityIndexes[2]] = 2;
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField& bfSf(mesh.Sf().boundaryField()[patchi]);
            vectorField bfNorm(bfSf/mag(bfSf));
            tmp<scalarField> scale;

            if (fixedWalls_[patchi])
            {
                scalarField m0(max(moments_(0).boundaryField()[patchi], 1e-8));
                tmp<scalarField> u(moments_(order100).boundaryField()[patchi]/m0);
                vectorField T(bfSf.size(), Zero);
                T.replace
                (
                    0,
                    max
                    (
                        moments_(order200).boundaryField()[patchi]/m0
                      - sqr(u),
                        1e-8
                    )
                );
                if (nd > 1)
                {
                    tmp<scalarField> v
                    (
                        moments_(order010).boundaryField()[patchi]/m0
                    );
                    T.replace
                    (
                        1,
                        max
                        (
                            moments_(order020).boundaryField()[patchi]/m0
                          - sqr(v),
                            1e-8
                        )
                    );
                }
                if (nd > 2)
                {
                    tmp<scalarField> w
                    (
                        moments_(order001).boundaryField()[patchi]/m0
                    );
                    T.replace
                    (
                        2,
                        max
                        (
                            moments_(order002).boundaryField()[patchi]/m0
                          - sqr(w),
                            1e-8
                        )
                    );
                }
                scale =
                (
                    sqrt
                    (
                        wallTemperatures_[patchi]*scalar(nd)
                       /(T & vector(1.0, 1.0, 1.0))
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

                if (scale.valid())
                {
                    bfUNei *= scale();
                }
                Gin += max(0.0, bfUOwn & bfSf)*bfwOwn;
                Gout -= min(0.0, bfUNei & bfSf)*bfwNei;
            }

            scalarField weightScale(Gin/(Gout + small));

            forAll(nodes, nodei)
            {
                scalarField& bfWNei =
                    nodesNei[nodei].primaryWeight().boundaryFieldRef()[patchi];
                bfWNei *= weightScale;
            }
        }
        else if (isA<symmetryFvPatch>(currPatch))
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

                bfwOwn = weight.boundaryField()[patchi].patchInternalField();
                bfwNei = bfwOwn;

                bfUOwn = U.boundaryField()[patchi].patchInternalField();
                bfUNei = bfUOwn - 2.0*(bfUOwn & bfNorm)*bfNorm;
            }
        }
    }
}


// ************************************************************************* //
