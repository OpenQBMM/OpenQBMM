/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Matteo Icardi
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

#include "growthModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(growthModel, 0);

    defineRunTimeSelectionTable(growthModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModel::growthModel
(
    const dictionary& dict
)
:
    dict_(dict),
    Cg_
    (
        dict.lookupOrDefault
        (
            "Cg",
            dimensionedScalar("one", inv(dimTime), 1.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModel::~growthModel()
{}


Foam::scalar
Foam::populationBalanceSubModels::growthModel::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    scalar gSource = 0.0;

    const PtrList<volNode>& nodes = quadrature.nodes();
    bool massBased = nodes[0].massBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return gSource;
    }
    label sizeOrder = momentOrder[sizeIndex];

    if (sizeOrder < 1)
    {
        return gSource;
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volNode& node = nodes[pNodeI];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], 0.0);

            scalar gSourcei = 2.0;
//             if (massBased)
//             {
//                 gSourcei = nodeSource(bAbscissa, sizeOrder);
//             }
//             else
//             {
//                 gSourcei = massNodeSource(bAbscissa, sizeOrder);
//             }

            gSourcei *=
                node.primaryWeight()[celli]
               *Kg(node.primaryAbscissae()[0][celli])
               *sizeOrder*pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    gSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }
            gSource += gSourcei;
        }

        return gSource;
    }

    forAll(nodes, pNodeI)
    {
        const volNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights()[0], sNodei)
        {
            scalar bAbscissa =
                max(node.secondaryAbscissae()[sizeIndex][sNodei][celli], 0.0);

//             scalar gSourcei = 0.0;
//             if (massBased)
//             {
//                 gSourcei = nodeSource(bAbscissa, sizeOrder);
//             }
//             else
//             {
//                 gSourcei = massNodeSource(bAbscissa, sizeOrder);
//             }

            scalar gSourcei =
                node.primaryWeight()[celli]
               *node.secondaryWeights()[sizeIndex][sNodei][celli]
               *Kg(bAbscissa)
               *sizeOrder*pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    gSourcei *=
                        node.secondaryWeights()[cmpt][sNodei][celli]
                       *pow
                        (
                            node.secondaryAbscissae()[cmpt][sNodei][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            gSource += gSourcei;
        }
    }

    return gSource;
}

Foam::scalar
Foam::populationBalanceSubModels::growthModel::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    scalar gSource = 0.0;

    const PtrList<volVelocityNode>& nodes = quadrature.nodes();
    bool massBased = nodes[0].massBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return gSource;
    }
    label sizeOrder = momentOrder[sizeIndex];

    if (sizeOrder < 1)
    {
        return gSource;
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volVelocityNode& node = nodes[pNodeI];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], 0.0);

//             scalar gSourcei = 0.0;
//             if (massBased)
//             {
//                 gSourcei = nodeSource(bAbscissa, sizeOrder);
//             }
//             else
//             {
//                 gSourcei = massNodeSource(bAbscissa, sizeOrder);
//             }

            scalar gSourcei =
                node.primaryWeight()[celli]
               *Kg(node.primaryAbscissae()[0][celli])
               *sizeOrder*pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    gSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }
            forAll(velocityIndexes, cmpt)
            {
                gSourcei *=
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
            gSource += gSourcei;
        }

        return gSource;
    }

    forAll(nodes, pNodeI)
    {
        const volVelocityNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights()[0], sNodei)
        {
            scalar bAbscissa =
                max(node.secondaryAbscissae()[sizeIndex][sNodei][celli], 0.0);

//             scalar gSourcei = 0.0;
//             if (massBased)
//             {
//                 gSourcei = nodeSource(bAbscissa, sizeOrder);
//             }
//             else
//             {
//                 gSourcei = massNodeSource(bAbscissa, sizeOrder);
//             }

            scalar gSourcei =
                node.primaryWeight()[celli]
               *node.secondaryWeights()[sizeIndex][sNodei][celli]
               *Kg(bAbscissa)
               *sizeOrder*pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    gSourcei *=
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
                gSourcei *=
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
            gSource += gSourcei;
        }
    }

    return gSource;
}

// ************************************************************************* //
