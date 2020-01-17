/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
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
    bool lengthBased = nodes[0].lengthBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return aSource;
    }

    label sizeOrder = momentOrder[sizeIndex];
    bool volumeFraction = nodes[0].useVolumeFraction();
    if (volumeFraction)
    {
        if (lengthBased)
        {
            sizeOrder += 3;
        }
        else
        {
            sizeOrder += 1;
        }
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();

    if (!nodes[0].extended())   // Non-extended quadrature case
    {
        forAll(nodes, pNode1i)
        {
            const volScalarNode& node1 = nodes[pNode1i];
            const volScalarField& pWeight1 = node1.primaryWeight();
            const PtrList<volScalarField>& pAbscissae1 = node1.primaryAbscissae();
            scalar bAbscissa1 =
                max(node1.primaryAbscissae()[sizeIndex][celli], scalar(0));
            scalar d1 = node1.d(celli, bAbscissa1);
            scalar n1 = node1.n(celli, pWeight1[celli], bAbscissa1);

            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();
                // Remove SMALL negative values in abscissae
                scalar bAbscissa2 =
                    max(node2.primaryAbscissae()[sizeIndex][celli], scalar(0));
                scalar d2 = node2.d(celli, bAbscissa2);
                scalar n2 = node2.n(celli, pWeight2[celli], bAbscissa2);


                scalar aSrc = 0.0;
                if (lengthBased)
                {
                    aSrc = nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }
                else
                {
                    aSrc = massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }

                aSrc *= n1*n2*Ka(d1, d2, Zero, celli, enviroment);

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
            const scalar bAbscissa1 =
                max
                (
                    node1.secondaryAbscissae()[sizeIndex][sNode1i][celli],
                    scalar(0)
                );
            scalar d1 = node1.d(celli, bAbscissa1);
            scalar n1 =
                node1.n(celli, pWeight1[celli], bAbscissa1)
               *node1.secondaryWeights()[sizeIndex][sNode1i][celli];

            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();

                forAll(node2.secondaryWeights()[0], sNode2i)
                {
                    scalar bAbscissa2 =
                        max
                        (
                            node2.secondaryAbscissae()[sizeIndex][sNode2i][celli],
                            scalar(0)
                        );
                    scalar d2 = node2.d(celli, bAbscissa2);
                    scalar n2 =
                        node2.n(celli, pWeight2[celli], bAbscissa2)
                       *node2.secondaryWeights()[sizeIndex][sNode2i][celli];

                    scalar aSrc = 0.0;
                    if (lengthBased)
                    {
                        aSrc =
                            nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                    }
                    else
                    {
                        aSrc =
                            massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                    }

                    aSrc *= n1*n2*Ka(d1, d2, Zero, celli, enviroment);
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

    const mappedPtrList<volVelocityNode>& nodes = quadrature.nodes();
    bool lengthBased = nodes[0].lengthBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return aSource;
    }

    label sizeOrder = momentOrder[sizeIndex];
    bool volumeFraction = nodes[0].useVolumeFraction();
    if (volumeFraction)
    {
        if (lengthBased)
        {
            sizeOrder += 3;
        }
        else
        {
            sizeOrder += 1;
        }
    }

    label nSizes = quadrature.nNodes()[sizeIndex];
    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    if (!pureSize())   // Ka is not only a function of size
    {
        forAll(nodes, pNode1i)
        {
            const volVelocityNode& node1 = nodes[pNode1i];
            scalar pWeight1 = node1.primaryWeight()[celli];
            scalar bAbscissa1 =
                max(node1.primaryAbscissae()[sizeIndex][celli], scalar(0));
            scalar d1 = node1.d(celli, bAbscissa1);
            scalar n1 = node1.n(celli, pWeight1, bAbscissa1);
            vector U1 = node1.velocityAbscissae()[celli];

            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volVelocityNode& node2 = nodes[pNode2i];
                scalar pWeight2 = node2.primaryWeight()[celli];
                scalar bAbscissa2 =
                    max(node2.primaryAbscissae()[sizeIndex][celli], scalar(0));
                scalar d2 = node2.d(celli, bAbscissa2);
                scalar n2 = node2.n(celli, pWeight2, bAbscissa2);
                vector U2 = node2.velocityAbscissae()[celli];

                scalar aSrc = 0.0;
                if (lengthBased)
                {
                    aSrc = nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }
                else
                {
                    aSrc = massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }

                aSrc *= n1*n2*Ka(d1, d2, U1 - U2, celli, enviroment);
                aSourcei += aSrc;

            }
            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    aSourcei *=
                        pow
                        (
                            node1.primaryAbscissae()[cmpt][celli],
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

    if (!nodes[0].extended())   // Non-extended quadrature case
    {
        scalarListList aSources(nSizes, scalarList(nSizes, Zero));
        for (label sizei = 0; sizei < nSizes; sizei++)
        {
            const volVelocityNode& node1 = nodes(sizei);
            scalar bAbscissa1 =
                max(node1.primaryAbscissae()[sizeIndex][celli], scalar(0));
            scalar d1 = node1.d(celli, bAbscissa1);

            for (label sizej = 0; sizej < nSizes; sizej++)
            {
                const volVelocityNode& node2 = nodes(sizej);
                scalar bAbscissa2 =
                    max(node2.primaryAbscissae()[sizeIndex][celli], scalar(0));
                scalar d2 = node2.d(celli, bAbscissa2);

                if (lengthBased)
                {
                    aSources[sizei][sizej] =
                        nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }
                else
                {
                    aSources[sizei][sizej] =
                        massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                }

                aSources[sizei][sizej] *=
                    Ka(d1, d2, Zero, celli, enviroment);

                if (volumeFraction)
                {
                    if (lengthBased)
                    {
                        aSources[sizei][sizej] /=
                            pow3(max(bAbscissa1, SMALL))
                           *pow3(max(bAbscissa2, SMALL));
                    }
                    else
                    {
                        aSources[sizei][sizej] /=
                            max(bAbscissa1, SMALL)
                           *max(bAbscissa2, SMALL) ;
                    }
                }
            }
        }

        forAll(nodes, pNode1i)
        {
            const volVelocityNode& node1 = nodes[pNode1i];
            label sizei = quadrature.nodeIndexes()[pNode1i][sizeIndex];
            scalar aSourcei = 0.0;

            forAll(nodes, pNode2i)
            {
                const volVelocityNode& node2 = nodes[pNode2i];
                label sizej = quadrature.nodeIndexes()[pNode2i][sizeIndex];

                aSourcei +=
                    node1.primaryWeight()[celli]
                   *node2.primaryWeight()[celli]
                   *aSources[sizei][sizej];

            }
            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    aSourcei *=
                        pow
                        (
                            node1.primaryAbscissae()[cmpt][celli],
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

    // Extended quadrature case
    label nSecondaryNodes(nodes[0].secondaryWeights().size());
    List<List<scalarListList>> aSources
    (
        nSizes,
        List<scalarListList>
        (
            nSecondaryNodes,
            scalarListList(nSizes, scalarList(nSecondaryNodes, Zero))
        )
    );
    for (label pNode1i = 0; pNode1i < nSizes; pNode1i++)
    {
        const volVelocityNode& node1 = nodes(pNode1i);
        for (label sNode1i = 0; sNode1i < nSecondaryNodes; sNode1i++)
        {
            scalar bAbscissa1 =
                max
                (
                    node1.secondaryAbscissae()[sizeIndex][sNode1i][celli],
                    scalar(0)
                );
            scalar d1 = node1.d(celli, bAbscissa1);

            for (label pNode2i = 0; pNode2i < nSizes; pNode2i++)
            {
                const volVelocityNode& node2 = nodes(pNode2i);
                for (label sNode2i = 0; sNode2i < nSecondaryNodes; sNode2i++)
                {
                    scalar bAbscissa2 =
                        max
                        (
                            node2.secondaryAbscissae()[sizeIndex][sNode2i][celli],
                            scalar(0)
                        );
                    scalar d2 = node2.d(celli, bAbscissa2);

                    if (lengthBased)
                    {
                        aSources[pNode1i][sNode1i][pNode2i][sNode2i] =
                            nodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                    }
                    else
                    {
                        aSources[pNode1i][sNode1i][pNode2i][sNode2i] =
                            massNodeSource(bAbscissa1, bAbscissa2, sizeOrder);
                    }

                    aSources[pNode1i][sNode1i][pNode2i][sNode2i] *=
                        Ka(d1, d2, Zero, celli, enviroment);

                    if (volumeFraction)
                    {
                        if (lengthBased)
                        {
                            aSources[pNode1i][sNode1i][pNode2i][sNode2i] /=
                                pow3(max(bAbscissa1, SMALL))
                               *pow3(max(bAbscissa2, SMALL));
                        }
                        else
                        {
                            aSources[pNode1i][sNode1i][pNode2i][sNode2i] /=
                                max(bAbscissa1, SMALL)
                               *max(bAbscissa2, SMALL);
                        }
                    }
                }
            }
        }
    }

    forAll(nodes, pNode1i)
    {
        const volVelocityNode& node1 = nodes[pNode1i];
        const volScalarField& pWeight1 = node1.primaryWeight();
        label sizei = quadrature.nodeIndexes()[pNode1i][sizeIndex];

        scalar aSourcei = 0.0;
        forAll(node1.secondaryWeights()[0], sNode1i)
        {
            forAll(nodes, pNode2i)
            {
                const volVelocityNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();
                label sizej = quadrature.nodeIndexes()[pNode2i][sizeIndex];

                forAll(node2.secondaryWeights()[0], sNode2i)
                {
                    aSourcei +=
                        pWeight1[celli]
                       *node1.secondaryWeights()[sizeIndex][sNode1i][celli]
                       *pWeight2[celli]
                       *node2.secondaryWeights()[sizeIndex][sNode2i][celli]
                       *aSources[sizei][sNode1i][sizej][sNode2i];
                }
            }
        }
        forAll(scalarIndexes, cmpt)
        {
            if (scalarIndexes[cmpt] != sizeIndex)
            {
                aSourcei *=
                    pow
                    (
                        node1.primaryAbscissae()[cmpt][celli],
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

// ************************************************************************* //
