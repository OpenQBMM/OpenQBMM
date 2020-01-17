/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
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

Foam::scalar
Foam::populationBalanceSubModels::collisionKernel::d
(
    const label nodei,
    const label celli
) const
{
    if (sizeIndex_ == -1)
    {
        return dp_()[celli];
    }

    const volVelocityNode& node = quadrature_.nodes()(nodei);
    scalar abscissa = node.primaryAbscissae()[sizeIndex_][celli];

    if (node.lengthBased())
    {
        return max(abscissa, minD_);
    }
    else
    {
        return cbrt(abscissa/rhos_[nodei]/Foam::constant::mathematical::pi*6.0);
    }
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
    nSizes_(0),
    rhos_(1),
    minD_(dict.lookupOrDefault("minD", SMALL)),
    implicit_(dict_.lookupOrDefault("implicit", true))
{
    if (sizeIndex_ != -1)
    {
        forAll(nodeIndexes_, nodei)
        {
            nSizes_ = max(nSizes_, nodeIndexes_[nodei][sizeIndex_] + 1);
        }

        rhos_.resize(nSizes_);
    }
    else
    {
        dp_=
            lookupOrInitialize
            (
                mesh,
                IOobject::groupName("d", quadrature.moments()[0].group()),
                dict,
                "d",
                dimLength
            ).ptr();
    }

    scalarList rhos(dict_.lookupOrDefault("rhos", scalarList()));

    forAll(rhos, i)
    {
        rhos_[i] = rhos[i];
    }

    if (rhos.size() < nSizes_ || nSizes_ == 0)
    {
        tmp<volScalarField> rho
        (
            lookupOrInitialize
            (
                mesh,
                IOobject::groupName
                (
                    "thermo:rho",
                    quadrature.moments()[0].group()
                ),
                dict,
                "rho",
                dimDensity
            )
        );

        for (label i = rhos.size(); i < nSizes_; i++)
        {
            rhos_[i] = rho()[0];
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

    if (nodeIndexes_[0].size() >= velocityIndexes_.size())
    {
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
