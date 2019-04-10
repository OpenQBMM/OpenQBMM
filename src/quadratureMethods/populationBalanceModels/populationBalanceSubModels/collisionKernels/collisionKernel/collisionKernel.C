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

#include "collisionKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(collisionKernel, 0);

    defineRunTimeSelectionTable(collisionKernel, dictionary);
}
}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::collisionKernel::lookupOrInitialize
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& dict,
    const word& entryName,
    const dimensionSet& dims
)
{
    if (mesh.foundObject<volScalarField>(name))
    {
        return mesh.lookupObject<volScalarField>(name);
    }
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(entryName, dims, dict)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernel::collisionKernel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature
)
:
    dict_(dict),
    mesh_(mesh),
    quadrature_(quadrature),
    momentOrders_(quadrature.momentOrders()),
    nodeIndexes_(quadrature.nodeIndexes()),
    velocityIndexes_(quadrature.nodes()[0].velocityIndexes()),
    nDimensions_(velocityIndexes_.size()),
    sizeIndex_(quadrature.nodes()[0].sizeIndex()),
    implicit_(dict_.lookupOrDefault("implicit", true))
{
    if (sizeIndex_ != -1)
    {
        forAll(nodeIndexes_, nodei)
        {
            nSizes_ = max(nSizes_, nodeIndexes_[nodei][sizeIndex_] + 1);
        }
    }

    forAll(momentOrders_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        labelList order(nDimensions_, 0);
        forAll(velocityIndexes_, cmpt)
        {
            order[cmpt] = momentOrder[velocityIndexes_[cmpt]];
        }

        bool found = false;
        forAll(velocityMomentOrders_, vmi)
        {
            bool same = true;
            forAll(velocityMomentOrders_[vmi], cmpt)
            {
                if (velocityMomentOrders_[vmi][cmpt] != order[cmpt])
                {
                    same = false;
                }
            }
            if (same)
            {
                found = true;
            }
        }
        if (!found)
        {
            velocityMomentOrders_.append(order);
        }
    }

    forAll(nodeIndexes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        labelList index(nDimensions_, 0);
        forAll(velocityIndexes_, cmpt)
        {
            index[cmpt] = nodeIndex[velocityIndexes_[cmpt]];
        }

        bool found = false;
        forAll(velocityNodeIndexes_, vmi)
        {
            bool same = true;
            forAll(velocityNodeIndexes_[vmi], cmpt)
            {
                if (velocityNodeIndexes_[vmi][cmpt] != index[cmpt])
                {
                    same = false;
                }
            }
            if (same)
            {
                found = true;
            }
        }
        if (!found)
        {
            velocityNodeIndexes_.append(index);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernel::~collisionKernel()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernel::preUpdate()
{}


void Foam::populationBalanceSubModels::collisionKernel::updateFields()
{
    if (!implicit_)
    {
        return;
    }

    forAll(quadrature_.moments()[0], celli)
    {
        updateCells(celli);
    }
}


// ************************************************************************* //
