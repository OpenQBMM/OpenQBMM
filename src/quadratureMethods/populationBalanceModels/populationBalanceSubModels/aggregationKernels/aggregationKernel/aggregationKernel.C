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

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(aggregationKernel, 0);

    defineRunTimeSelectionTable(aggregationKernel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernel::aggregationKernel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    Ca_
    (
        dict.lookupOrDefault
        (
            "Ca",
            dimensionedScalar("one", inv(dimTime), 1.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernel::~aggregationKernel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernel::nodeSource
(
    const scalar& abscissa1,
    const scalar& abscissa2,
    const label momentOrder
) const
{
    return
        0.5
       *pow // Birth
        (
            pow3(abscissa1) + pow3(abscissa2),
            momentOrder/3.0
        )
      - pow(abscissa1, momentOrder);
}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernel::massNodeSource
(
    const scalar& abscissa1,
    const scalar& abscissa2,
    const label momentOrder
) const
{
    return
        0.5
       *pow // Birth
        (
            abscissa1 + abscissa2,
            momentOrder
        )
      - pow(abscissa1, momentOrder);
}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernel::aggregationSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label enviroment
)
{
    scalar aSource = 0.0;

    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    bool massBased = nodes[0].massBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return aSource;
    }

    label sizeOrder = momentOrder[sizeIndex];

    const labelList& scalarIndexes = nodes[0].scalarIndexes();

    if (!nodes[0].extended())   // Non-extended quadrature case
    {
        forAll(nodes, pNode1i)
        {
            const volScalarNode& node1 = nodes[pNode1i];
            const volScalarField& pWeight1 = node1.primaryWeight();
            const PtrList<volScalarField>& pAbscissae1 =
                node1.primaryAbscissae();
            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();

                // Remove small negative values in abscissae
                scalar bAbscissa1 =
                    max
                    (
                        node1.primaryAbscissae()[sizeIndex][celli],
                        0.0
                    );
                scalar bAbscissa2 =
                    max
                    (
                        node2.primaryAbscissae()[sizeIndex][celli],
                        0.0
                    );

                scalar aSrc = 0.0;
                if (massBased)
                {
                    aSrc =
                        massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }
                else
                {
                    aSrc = nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }

                aSrc *=
                    pWeight1[celli]
                   *pWeight2[celli]
                   *Ka(bAbscissa1, bAbscissa2, celli, enviroment);

                aSourcei += aSrc;
            }

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    aSourcei *=
                        pow
                        (
                            pAbscissae1[cmpt][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            aSource += aSourcei;
        }

        return aSource;
    }

    forAll(nodes, pNode1i)      // Extended quadrature case
    {
        const volScalarNode& node1 = nodes[pNode1i];
        const volScalarField& pWeight1 = node1.primaryWeight();

        forAll(node1.secondaryWeights()[0], sNode1i)
        {
            const volScalarField& sAbscissa1 =
                node1.secondaryAbscissae()[sizeIndex][sNode1i];

            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();

                forAll(node2.secondaryWeights()[0], sNode2i)
                {
                    const volScalarField& sAbscissa2
                        = node2.secondaryAbscissae()[sizeIndex][sNode2i];

                    // Remove small negative values in abscissae
                    scalar bAbscissa1 = max(sAbscissa1[celli], 0.0);
                    scalar bAbscissa2 = max(sAbscissa2[celli], 0.0);

                    scalar aSrc = 0.0;
                    if (massBased)
                    {
                        aSrc =
                            massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                    }
                    else
                    {
                        aSrc = nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                    }

                    aSrc *=
                        pWeight1[celli]
                       *node1.secondaryWeights()[sizeIndex][sNode1i][celli]
                       *pWeight2[celli]
                       *node2.secondaryWeights()[sizeIndex][sNode2i][celli]
                       *Ka(bAbscissa1, bAbscissa2, celli, enviroment);

                    aSourcei += aSrc;
                }
            }

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    aSourcei *=
                        node1.secondaryWeights()[cmpt][sNode1i][celli]
                       *pow
                        (
                            node1.secondaryAbscissae()[cmpt][sNode1i][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            aSource += aSourcei;
        }
    }

    return aSource;
}

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernel::aggregationSource
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature,
    const label enviroment
)
{
    scalar aSource = 0.0;

    const PtrList<volVelocityNode>& nodes = quadrature.nodes();
    bool massBased = nodes[0].massBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return aSource;
    }

    label sizeOrder = momentOrder[sizeIndex];

    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    if (!nodes[0].extended())   // Non-extended quadrature case
    {
        forAll(nodes, pNode1i)
        {
            const volScalarNode& node1 = nodes[pNode1i];
            const volScalarField& pWeight1 = node1.primaryWeight();
            const PtrList<volScalarField>& pAbscissae1 =
                node1.primaryAbscissae();
            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();

                // Remove small negative values in abscissae
                scalar bAbscissa1 =
                    max
                    (
                        node1.primaryAbscissae()[sizeIndex][celli],
                        0.0
                    );
                scalar bAbscissa2 =
                    max
                    (
                        node2.primaryAbscissae()[sizeIndex][celli],
                        0.0
                    );

                scalar aSrc = 0.0;
                if (massBased)
                {
                    aSrc =
                        massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }
                else
                {
                    aSrc = nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }

                aSrc *=
                    pWeight1[celli]
                   *pWeight2[celli]
                   *Ka(bAbscissa1, bAbscissa2, celli, enviroment);

                aSourcei += aSrc;
            }
            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    aSourcei *=
                        pow
                        (
                            pAbscissae1[cmpt][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            forAll(velocityIndexes, cmpt)
            {
                aSourcei *=
                    pow
                    (
                        node1.velocityAbscissae()[celli][cmpt],
                        momentOrder[velocityIndexes[cmpt]]
                    );
            }
            aSource += aSourcei;
        }

        return aSource;
    }

    forAll(nodes, pNode1i)      // Extended quadrature case
    {
        const volVelocityNode& node1 = nodes[pNode1i];
        const volScalarField& pWeight1 = node1.primaryWeight();

        forAll(node1.secondaryWeights()[0], sNode1i)
        {
            const volScalarField& sAbscissa1 =
                node1.secondaryAbscissae()[sizeIndex][sNode1i];

            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();

                forAll(node2.secondaryWeights()[0], sNode2i)
                {
                    const volScalarField& sAbscissa2
                        = node2.secondaryAbscissae()[sizeIndex][sNode2i];

                    // Remove small negative values in abscissae
                    scalar bAbscissa1 = max(sAbscissa1[celli], 0.0);
                    scalar bAbscissa2 = max(sAbscissa2[celli], 0.0);

                    scalar aSrc = 0.0;
                    if (massBased)
                    {
                        aSrc =
                            massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                    }
                    else
                    {
                        aSrc =
                            nodeSource
                            (
                                bAbscissa1,
                                bAbscissa2,
                                sizeOrder
                            );
                    }

                    aSrc *=
                        pWeight1[celli]
                       *node1.secondaryWeights()[sizeIndex][sNode1i][celli]
                       *pWeight2[celli]
                       *node2.secondaryWeights()[sizeIndex][sNode2i][celli]
                       *Ka(bAbscissa1, bAbscissa2, celli, enviroment);

                    aSourcei += aSrc;
                }
            }
            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    aSourcei *=
                        node1.secondaryWeights()[cmpt][sNode1i][celli]
                       *pow
                        (
                            node1.secondaryAbscissae()[cmpt][sNode1i][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            forAll(velocityIndexes, cmpt)
            {
                aSourcei *=
                    pow
                    (
                        node1.velocityAbscissae()[celli][cmpt],
                        momentOrder[velocityIndexes[cmpt]]
                    );
            }
            aSource += aSourcei;
        }
    }

    return aSource;
}

// ************************************************************************* //
