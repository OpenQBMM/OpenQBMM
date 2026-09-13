/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2026 Alberto Passalacqua
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

#include "TensorProductMomentInversion.H"
#include "addToRunTimeSelectionTable.H"
#include "supportType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(TensorProduct, 0);
    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        TensorProduct,
        dictionary
    );
}
}

void Foam::multivariateMomentInversions::TensorProduct::buildIndexes
(
    labelListList& nodeIndexes,
    const labelList& nNodes,
    label dimi,
    label& nodei,
    labelList& index
)
{
    if (dimi < nNodes.size())
    {
        for (label i = 0; i < nNodes[dimi]; i++)
        {
            index[dimi] = i;
            buildIndexes(nodeIndexes, nNodes, dimi+1, nodei, index);
        }
    }
    else
    {
        nodeIndexes[nodei] = index;
        nodei++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::TensorProduct::TensorProduct
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
:
    multivariateMomentInversion
    (
        dict,
        momentOrders,
        nodeIndexes,
        velocityIndexes
    ),
    nPureMoments_(nNodes_.size(), 0),
    supports_(wordListToSupportTypeList(dict.lookup("supports"))),
    univariateInverters_(nNodes_.size()),
    smallM0_(VGREAT),
    smallZeta_(VGREAT)
{
    forAll(univariateInverters_, dimi)
    {
        univariateInverters_.set
        (
            dimi,
            univariateMomentInversion::New
            (
                dict.subDict("basicQuadrature" + Foam::name(dimi))
            ).ptr()
        );

        // Each direction is inverted on its own, by its own quadrature, and
        // with its own thresholds. What the inversion reports as its
        // thresholds is the smallest of them: the directions share m0, so
        // a check made before any of them is inverted should be no stricter
        // than the most permissive, and each still applies its own.
        //
        // They used to be the largest, floored at SMALL, and imposed on every
        // direction. A threshold could then only ever be raised, and the
        // zeta_k, which carry the units of the abscissa of their direction,
        // were cut at 1e-15 whatever those units were.
        smallM0_ = min(smallM0_, univariateInverters_[dimi].smallM0());
        smallZeta_ = min(smallZeta_, univariateInverters_[dimi].smallZeta());
    }

    forAll(momentOrders_, mi)
    {
        forAll(nPureMoments_, dimi)
        {
            nPureMoments_[dimi] =
                max
                (
                    nPureMoments_[dimi],
                    momentOrders_[mi][dimi] + 1
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::TensorProduct::~TensorProduct()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multivariateMomentInversions::TensorProduct::invert
(
    const multivariateMomentSet& moments
)
{
    reset();
    labelList nNonZeroNodes(nNodes_.size(), 0);
    labelList zeroOrder(momentOrders_[0].size(), 0);

    // The largest abscissa of each direction, which the system for the
    // weights is solved in units of below
    scalarList directionScale(nNodes_.size(), 1.0);

    label vi = 0;
    label si = 0;

    forAll(univariateInverters_, dimi)
    {
        univariateMomentSet univariateMoments
        (
            nPureMoments_[dimi],
            supports_[dimi],
            univariateInverters_[dimi].smallM0(),
            univariateInverters_[dimi].smallZeta(),
            Zero
        );

        forAll(univariateMoments, mi)
        {
            labelList momentOrder(zeroOrder);
            momentOrder[dimi] = mi;
            univariateMoments[mi] = moments(momentOrder);
        }

        if (!univariateMoments.isRealizable(false))
        {
            return false;
        }

        univariateInverters_[dimi].invert(univariateMoments);

        const scalarList& abscissae =
            univariateInverters_[dimi].abscissae();

        nNonZeroNodes[dimi] = abscissae.size();

        scalar largestAbscissa = 0;

        forAll(abscissae, nodei)
        {
            largestAbscissa = max(largestAbscissa, mag(abscissae[nodei]));
        }

        // A direction whose nodes all sit at zero is left in its own units,
        // in which the system is not affected by it
        directionScale[dimi] = largestAbscissa > 0 ? largestAbscissa : 1.0;

        if
        (
            max(univariateInverters_[dimi].weights())
          < univariateInverters_[dimi].smallM0()
        )
        {
            nNonZeroNodes[dimi] = 0;
        }

        if (nNonZeroNodes[dimi] > 0)
        {
            forAll(nodeIndexes_, nodei)
            {
                label nodeIndex = nodeIndexes_[nodei][dimi];

                if (nodeIndex < abscissae.size())
                {
                    if (isVelocityDimension(dimi, vi))
                    {
                        velocityAbscissae_[nodei][vi] = abscissae[nodeIndex];
                    }
                    else
                    {
                        abscissae_[nodei][si] = abscissae[nodeIndex];
                    }
                }
            }
        }

        if (isVelocityDimension(dimi, vi))
        {
            vi++;
        }
        else
        {
            si++;
        }
    }

    // The quadrature has a node for every combination of the univariate
    // abscissae, so a direction that carries no node leaves it with none at
    // all, whatever the other directions carry. The weights are left at the
    // zero reset() wrote.
    //
    // The test is on the smallest of the counts rather than the largest:
    // taking the product over the directions that do carry nodes, and
    // building the indexes over all of them, would ask buildIndexes for
    // more nodes than it writes and leave the rest of the list at zero,
    // which is a set of identical indexes and a singular system to solve
    // for the weights.
    if (min(nNonZeroNodes) == 0)
    {
        return true;
    }

    label totNonZeroNodes = 1;

    forAll(nNonZeroNodes, dimi)
    {
        totNonZeroNodes *= nNonZeroNodes[dimi];
    }

    const label nDims = nNonZeroNodes.size();

    labelListList nonZeroNodeIndexes(totNonZeroNodes, labelList(nDims, 0));
    {
        label nodei = 0;
        labelList index(nDims, 0);
        buildIndexes(nonZeroNodeIndexes, nNonZeroNodes, 0, nodei, index);
    }

    scalarList mixedMoments(nonZeroNodeIndexes.size(), Zero);
    scalarSquareMatrix R(nonZeroNodeIndexes.size(), 1.0);

    // The system is a tensor-product Vandermonde matrix, whose row of
    // orders (k, l, m) holds the abscissae of each node raised to those
    // orders. In the units of the case its entries span one to the scale of
    // a direction to the power of the largest order, and a direction small
    // enough in its own units, as a volume in cubic metres is, leaves it
    // singular to the elimination. Each abscissa is therefore divided by the
    // scale of its direction, and each mixed moment by the scales raised to
    // its orders, which divides both sides of the system by the same factor
    // row by row and leaves the weights unchanged.
    forAll(nonZeroNodeIndexes, nodei)
    {
        const labelList& order = nonZeroNodeIndexes[nodei];

        scalar scaleOfOrder = 1.0;

        forAll(order, dimi)
        {
            scaleOfOrder *= pow(directionScale[dimi], order[dimi]);
        }

        mixedMoments[nodei] = moments(order)/scaleOfOrder;
    }

    forAll(nonZeroNodeIndexes, mi)
    {
        forAll(nonZeroNodeIndexes, nodei)
        {
            vi = 0;
            si = 0;

            forAll(nonZeroNodeIndexes[nodei], dimi)
            {
                if (isVelocityDimension(dimi, vi))
                {
                    R(mi, nodei) *=
                        pow
                        (
                            velocityAbscissae_(nonZeroNodeIndexes[nodei])[vi]
                           /directionScale[dimi],
                            nonZeroNodeIndexes[mi][dimi]
                        );
                    vi++;
                }
                else
                {
                    R(mi, nodei) *=
                        pow
                        (
                            abscissae_(nonZeroNodeIndexes[nodei])[si]
                           /directionScale[dimi],
                            nonZeroNodeIndexes[mi][dimi]
                        );
                    si++;
                }
            }
        }
    }

    // One weight per node, not per direction
    scalarList weights(totNonZeroNodes, Zero);
    solve(weights, R, mixedMoments);

    forAll(nonZeroNodeIndexes, nodei)
    {
        weights_(nonZeroNodeIndexes[nodei]) = weights[nodei];
    }

    return true;
}

Foam::scalar Foam::multivariateMomentInversions::TensorProduct::smallM0() const
{
    return smallM0_;
}

Foam::scalar Foam::multivariateMomentInversions::TensorProduct::smallZeta() const
{
    return smallZeta_;
}

// ************************************************************************* //
