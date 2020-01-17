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
            dict.subDict("daughterDistribution")
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
    bool lengthBased = nodes[0].lengthBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return bSource;
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


    if (!nodes[0].extended())
    {
        forAll(nodes, pNodei)
        {
             const volScalarNode& node = nodes[pNodei];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));

            scalar bSourcei = 0.0;

            if (lengthBased)
            {
                bSourcei = nodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSourcei = massNodeSource(bAbscissa, sizeOrder);
            }

            bSourcei *=
                node.primaryWeight()[celli]
               *Kb(bAbscissa, celli);

            if (volumeFraction)
            {
                if (lengthBased)
                {
                    bSourcei /= pow3(max(bAbscissa, SMALL));
                }
                else
                {
                    bSourcei /= max(bAbscissa, SMALL);
                }
            }

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
                max
                (
                    node.secondaryAbscissae()[sizeIndex][sNodei][celli],
                    scalar(0)
                );

            scalar bSourcei = 0.0;

            if (lengthBased)
            {
                bSourcei = nodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSourcei = massNodeSource(bAbscissa, sizeOrder);
            }

            bSourcei *=
                node.primaryWeight()[celli]
               *node.secondaryWeights()[sizeIndex][sNodei][celli]
               *Kb(bAbscissa, celli);

            if (volumeFraction)
            {
                if (lengthBased)
                {
                    bSourcei /= pow3(max(bAbscissa, SMALL));
                }
                else
                {
                    bSourcei /= max(bAbscissa, SMALL);
                }
            }

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

    const mappedPtrList<volVelocityNode>& nodes = quadrature.nodes();
    bool lengthBased = nodes[0].lengthBased();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return bSource;
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

    if (!nodes[0].extended())
    {
        scalarList bSources(nSizes, Zero);

        for (label sizei = 0; sizei < nSizes; sizei++)
        {
            const volVelocityNode& node = nodes(sizei);

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));

            if (lengthBased)
            {
                bSources[sizei] = nodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSources[sizei] = massNodeSource(bAbscissa, sizeOrder);
            }

            bSources[sizei] *= Kb(bAbscissa, celli);

            if (volumeFraction)
            {
                if (lengthBased)
                {
                    bSources[sizei] /= pow3(max(bAbscissa, SMALL));
                }
                else
                {
                    bSources[sizei] /= max(bAbscissa, SMALL);
                }
            }
        }

        forAll(nodes, pNodei)
        {
             const volVelocityNode& node = nodes[pNodei];
             label sizei = quadrature.nodeIndexes()[pNodei][sizeIndex];

            scalar bSourcei = node.primaryWeight()[celli]*bSources[sizei];

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
                        node.velocityAbscissae()[celli][cmpt],
                        momentOrder[velocityIndexes[cmpt]]
                    );
            }

            bSource += bSourcei;
        }

        return bSource;
    }

    label nSecondaryNodes(nodes[0].secondaryWeights().size());
    scalarListList bSources(nSizes, scalarList(nSecondaryNodes, Zero));
    for (label pNodei = 0; pNodei < nSizes; pNodei++)
    {
        for (label sNodei = 0; sNodei < nSecondaryNodes; sNodei++)
        {
            const volVelocityNode& node = nodes(pNodei);

            scalar bAbscissa =
                max
                (
                    node.secondaryAbscissae()[sizeIndex][sNodei][celli],
                    scalar(0)
                );

            if (lengthBased)
            {
                bSources[pNodei][sNodei] = nodeSource(bAbscissa, sizeOrder);
            }
            else
            {
                bSources[pNodei][sNodei] = massNodeSource(bAbscissa, sizeOrder);
            }

            bSources[pNodei][sNodei] *= Kb(bAbscissa, celli);

            if (volumeFraction)
            {
                if (lengthBased)
                {
                    bSources[pNodei][sNodei] /= pow3(max(bAbscissa, SMALL));
                }
                else
                {
                    bSources[pNodei][sNodei] /= max(bAbscissa, SMALL);
                }
            }
        }
    }

    forAll(nodes, pNodei)
    {
        const volVelocityNode& node = nodes[pNodei];
        label sizei = quadrature.nodeIndexes()[pNodei][sizeIndex];

        scalar bSourcei = 0.0;

        forAll(node.secondaryWeights()[0], sNodei)
        {
            bSourcei +=
                node.primaryWeight()[celli]
               *node.secondaryWeights()[sizeIndex][sNodei][celli]
               *bSources[sizei][sNodei];

        }

        forAll(scalarIndexes, cmpt)
        {
            if (scalarIndexes[cmpt] != sizeIndex)
            {
                bSourcei *=
                    pow
                    (
                        node.primaryAbscissae()[cmpt][celli],
                        momentOrder[scalarIndexes[cmpt]]
                    );
            }
        }

        forAll(velocityIndexes, cmpt)
        {
            bSourcei *=
                pow
                (
                    node.velocityAbscissae()[celli][cmpt],
                    momentOrder[velocityIndexes[cmpt]]
                );
        }
        
        bSource += bSourcei;
    }

    return bSource;
}

// ************************************************************************* //
