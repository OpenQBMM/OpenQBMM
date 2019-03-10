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

#include "breakupKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(breakupKernel, 0);

    defineRunTimeSelectionTable(breakupKernel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernel::breakupKernel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    Cb_
    (
        dict.lookupOrDefault
        (
            "Cb",
            dimensionedScalar("one", inv(dimTime), 1.0)
        )
    ),
    daughterDistribution_
    (
        Foam::populationBalanceSubModels::daughterDistribution::New
        (
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernel::~breakupKernel()
{}


Foam::scalar
Foam::populationBalanceSubModels::breakupKernel::nodeSource
(
    const scalar& abscissa,
    const label momentOrder
) const
{
    return
         daughterDistribution_->mD(momentOrder, abscissa)//Birth
       - pow(abscissa, momentOrder);
}


Foam::scalar
Foam::populationBalanceSubModels::breakupKernel::massNodeSource
(
    const scalar& abscissa,
    const label momentOrder
) const
{
    return
        daughterDistribution_->mDMass(momentOrder, abscissa)//Birth
      - pow(abscissa, momentOrder);
}


Foam::scalar
Foam::populationBalanceSubModels::breakupKernel::breakupSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    scalar bSource = 0.0;

    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    bool massBased = nodes[0].massBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return bSource;
    }

    label sizeOrder = momentOrder[sizeIndex];

    const labelList& scalarIndexes = nodes[0].scalarIndexes();


    if (!nodes[0].extended())
    {
        forAll(nodes, pNodei)
        {
             const volScalarNode& node = nodes[pNodei];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], 0.0);

            scalar bSourcei = 0.0;
            if (massBased)
            {
                bSourcei = massNodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSourcei = nodeSource(bAbscissa, sizeOrder);
            }

            bSourcei *=
                node.primaryWeight()[celli]
                *Kb(bAbscissa, celli);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    bSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }
            bSource += bSourcei;
        }

        return bSource;
    }

    forAll(nodes, pNodei)
    {
        const volScalarNode& node = nodes[pNodei];

        forAll(node.secondaryWeights()[0], sNodei)
        {
            scalar bAbscissa =
                max(node.secondaryAbscissae()[sizeIndex][sNodei][celli], 0.0);

            scalar bSourcei = 0.0;
            if (massBased)
            {
                bSourcei = massNodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSourcei = nodeSource(bAbscissa, sizeOrder);
            }

            bSourcei *=
                node.primaryWeight()[celli]
               *node.secondaryWeights()[sizeIndex][sNodei][celli]
               *Kb(bAbscissa, celli);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    bSourcei *=
                        node.secondaryWeights()[cmpt][sNodei][celli]
                       *pow
                        (
                            node.secondaryAbscissae()[cmpt][sNodei][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            bSource += bSourcei;
        }
    }

    return bSource;
}

Foam::scalar
Foam::populationBalanceSubModels::breakupKernel::breakupSource
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    scalar bSource = 0.0;

    const PtrList<volVelocityNode>& nodes = quadrature.nodes();
    bool massBased = nodes[0].massBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return bSource;
    }

    label sizeOrder = momentOrder[sizeIndex];

    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodei)
        {
             const volVelocityNode& node = nodes[pNodei];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], 0.0);

            scalar bSourcei = 0.0;
            if (massBased)
            {
                bSourcei = massNodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSourcei = nodeSource(bAbscissa, sizeOrder);
            }

            bSourcei *=
                node.primaryWeight()[celli]
                *Kb(bAbscissa, celli);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    bSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }
            forAll(velocityIndexes, cmpt)
            {
                bSourcei *=
                    pow
                    (
                        component
                        (
                            node.velocityAbscissae()[celli],
                            cmpt
                        ),
                        momentOrder[velocityIndexes[cmpt]]
                    );
            }
            bSource += bSourcei;
        }

        return bSource;
    }

    forAll(nodes, pNodei)
    {
        const volVelocityNode& node = nodes[pNodei];

        forAll(node.secondaryWeights()[0], sNodei)
        {
            scalar bAbscissa =
                max(node.secondaryAbscissae()[sizeIndex][sNodei][celli], 0.0);

            scalar bSourcei = 0.0;
            if (massBased)
            {
                bSourcei = massNodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSourcei = nodeSource(bAbscissa, sizeOrder);
            }

            bSourcei *=
                node.primaryWeight()[celli]
               *node.secondaryWeights()[sizeIndex][sNodei][celli]
               *Kb(bAbscissa, celli);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    bSourcei *=
                        node.secondaryWeights()[cmpt][sNodei][celli]
                       *pow
                        (
                            node.secondaryAbscissae()[cmpt][sNodei][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            forAll(velocityIndexes, cmpt)
            {
                bSourcei *=
                    pow
                    (
                        component
                        (
                            node.velocityAbscissae()[celli],
                            cmpt
                        ),
                        momentOrder[velocityIndexes[cmpt]]
                    );
            }
            bSource += bSourcei;
        }
    }

    return bSource;
}

// ************************************************************************* //
