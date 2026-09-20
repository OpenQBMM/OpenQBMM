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

#include "conditionalMeansMomentInversion.H"
#include "mappedLists.H"
#include "supportType.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(conditionalMeans, 0);

    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        conditionalMeans,
        dictionary
    );

    // The name the velocity closure has always been selected by
    addNamedToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        conditionalMeans,
        dictionary,
        monoKinetic
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::conditionalMeans::conditionalMeans
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
    primaryIndex_
    (
        calcPrimaryIndex(momentOrders[0].size(), velocityIndexes_)
    ),
    scalarIndexes_(),
    scalarComponents_(),
    nPrimaryMoments_(calcNPrimaryMoments(momentOrders, primaryIndex_)),
    primaryInverter_
    (
        univariateMomentInversion::New
        (
            dict.subDict("basicQuadrature"),
            nodeIndexes.size()
        )
    ),
    nAdditionalQuadraturePoints_
    (
        primaryInverter_().nAdditionalQuadraturePoints(nPrimaryMoments_)
    ),
    nClipped_(0),
    clipToSupports_(true)
{
    // Every node carries a conditional mean, found from as many moments
    // <x^k phi> as there are nodes. The moments of a scalar that is linear
    // in the primary coordinate, which is what an inert scalar and the
    // mixing of two streams give, involve the moments of the primary
    // coordinate up to the order of the number of nodes. Beyond the ones
    // that are transported they are what the closure of the primary
    // inversion says, not what the distribution has, and the conditional
    // solve amplifies the difference.
    if (nodeIndexes.size() > max(nPrimaryMoments_ - 1, label(1)))
    {
        WarningInFunction
            << "The primary coordinate has " << nodeIndexes.size()
            << " nodes and " << nPrimaryMoments_ << " moments." << nl
            << "    The conditional means of a scalar linear in the primary "
            << "coordinate are exact" << nl
            << "    with up to " << max(nPrimaryMoments_ - 1, label(1))
            << " nodes: with N regular nodes, up to N - 1 additional ones."
            << nl
            << "    Beyond that they depend on moments of the primary "
            << "coordinate that are" << nl
            << "    not transported, which the closure of the quadrature "
            << "supplies." << nl
            << endl;
    }

    // The abscissae hold the dimensions that are not velocities, in their
    // order: the primary one first, since it is the first of them, and the
    // conditioned scalars after it
    label component = 0;

    for (label dimi = 0; dimi < nDistributionDims_; dimi++)
    {
        if (velocityIndexes_.found(dimi))
        {
            continue;
        }

        if (dimi != primaryIndex_)
        {
            scalarIndexes_.append(dimi);
            scalarComponents_.append(component);
        }

        component++;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::conditionalMeans::~conditionalMeans()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label
Foam::multivariateMomentInversions::conditionalMeans::calcPrimaryIndex
(
    const label nDistributionDims,
    const labelList& velocityIndexes
)
{
    boolList isVelocity(nDistributionDims, false);

    forAll(velocityIndexes, vi)
    {
        const label dimi = velocityIndexes[vi];

        // The callers pass {-1} for a distribution with no velocity
        if (dimi < 0)
        {
            continue;
        }

        if (dimi >= nDistributionDims)
        {
            FatalErrorInFunction
                << "The velocity index " << dimi << " is not a dimension of "
                << "the distribution." << nl
                << "    Number of dimensions: " << nDistributionDims << nl
                << exit(FatalError);
        }

        isVelocity[dimi] = true;
    }

    forAll(isVelocity, dimi)
    {
        if (!isVelocity[dimi])
        {
            return dimi;
        }
    }

    FatalErrorInFunction
        << "Every dimension of the distribution is a velocity, so there is "
        << "no coordinate to condition the others on." << nl
        << "    Velocity indexes: " << velocityIndexes << nl
        << exit(FatalError);

    return -1;
}


Foam::label
Foam::multivariateMomentInversions::conditionalMeans::calcNPrimaryMoments
(
    const labelListList& momentOrders,
    const label primaryIndex
)
{
    // Only the pure moments count: with more nodes than half of them, as
    // a generalized quadrature has, the conditional moments reach higher
    // powers of the primary coordinate than the pure ones do
    label maxOrder = 0;

    forAll(momentOrders, mi)
    {
        const labelList& order = momentOrders[mi];

        bool pure = true;

        forAll(order, dimi)
        {
            pure = pure && (dimi == primaryIndex || order[dimi] == 0);
        }

        if (pure)
        {
            maxOrder = max(maxOrder, order[primaryIndex]);
        }
    }

    return maxOrder + 1;
}


void
Foam::multivariateMomentInversions::conditionalMeans::momentOrdersAndNodes
(
    const label nPrimaryMoments,
    const label nNodes,
    const label nScalars,
    const label nVelocityDimensions,
    labelListList& momentOrders,
    labelListList& nodeIndexes
)
{
    const label nDims = 1 + nScalars + nVelocityDimensions;

    momentOrders.clear();
    nodeIndexes.clear();

    for (label k = 0; k < nPrimaryMoments; k++)
    {
        labelList order(nDims, Zero);
        order[0] = k;
        momentOrders.append(order);
    }

    for (label dimi = 1; dimi < nDims; dimi++)
    {
        for (label k = 0; k < nNodes; k++)
        {
            labelList order(nDims, Zero);
            order[0] = k;
            order[dimi] = 1;
            momentOrders.append(order);
        }
    }

    for (label nodei = 0; nodei < nNodes; nodei++)
    {
        labelList index(nDims, Zero);
        index[0] = nodei;
        nodeIndexes.append(index);
    }
}


bool Foam::multivariateMomentInversions::conditionalMeans::invert
(
    const multivariateMomentSet& moments
)
{
    reset();
    nClipped_ = 0;

    const List<supportType>& supports = moments.supports();

    univariateMomentSet primaryMoments
    (
        nPrimaryMoments_,
        supports[primaryIndex_],
        primaryInverter_().smallM0(),
        primaryInverter_().smallZeta(),
        Zero,
        nAdditionalQuadraturePoints_
    );

    forAll(primaryMoments, mi)
    {
        labelList primaryOrder(nDistributionDims_, 0);
        primaryOrder[primaryIndex_] = mi;

        primaryMoments[mi] = moments(primaryOrder);
    }

    if (!primaryMoments.isRealizable(false))
    {
        return false;
    }

    primaryInverter_->invert(primaryMoments);
    const scalarList& primaryWeights(primaryInverter_->weights());
    const scalarList& primaryAbscissae(primaryInverter_->abscissae());

    forAll(primaryWeights, nodei)
    {
        weights_[nodei] = primaryWeights[nodei];
        abscissae_[nodei][0] = primaryAbscissae[nodei];
    }

    const label nNodes = primaryWeights.size();

    if (nNodes == 0)
    {
        return true;
    }

    scalarSquareMatrix invVR(nNodes, Zero);

    // The conditional solve all the quadratures share. Every node that
    // carries anything stays in the system, so that the moments are shared
    // out exactly, and a node that carries less than smallM0 is solved for
    // but not read back: its row is zero. Only a node that carries nothing
    // is left out.
    labelList solvedNodes;
    scalarSquareMatrix invVRSolved;

    if
    (
        !Vandermonde::conditionalInverse
        (
            primaryAbscissae,
            primaryWeights,
            0,
            smallM0(),
            solvedNodes,
            invVRSolved
        )
    )
    {
        return false;
    }

    forAll(solvedNodes, i)
    {
        forAll(solvedNodes, k)
        {
            invVR(solvedNodes[i], k) = invVRSolved(i, k);
        }
    }

    // The conditional means of the coordinate of dimension dimi at every
    // node, from the moments <x^k phi>, k < nNodes
    auto conditionalMeansOf = [&](const label dimi)
    {
        labelList order(nDistributionDims_, 0);
        order[dimi] = 1;

        scalarRectangularMatrix M(nNodes, 1, 0);

        for (label nodei = 0; nodei < nNodes; nodei++)
        {
            order[primaryIndex_] = nodei;
            M(nodei, 0) = moments(order);
        }

        return scalarRectangularMatrix(invVR*M);
    };

    // The velocities. vi counts the velocity components, which is how they
    // are stored in velocityAbscissae_, while velocityIndexes_[vi] is the
    // dimension of the distribution that component belongs to, which is
    // how the moment orders are written.
    for (label vi = 0; vi < nVelocityDimensions_; vi++)
    {
        const scalarRectangularMatrix nu
        (
            conditionalMeansOf(velocityIndexes_[vi])
        );

        forAll(primaryWeights, nodei)
        {
            if (primaryWeights[nodei] > smallM0())
            {
                velocityAbscissae_[nodei][vi] = nu(nodei, 0);
            }
        }
    }

    // The scalars, each kept within its support
    forAll(scalarIndexes_, si)
    {
        const label dimi = scalarIndexes_[si];
        const supportType support = supports[dimi];

        const scalarRectangularMatrix phi(conditionalMeansOf(dimi));

        forAll(primaryWeights, nodei)
        {
            if (primaryWeights[nodei] <= smallM0())
            {
                continue;
            }

            scalar value = phi(nodei, 0);

            if (!clipToSupports_)
            {
                // The model the quadrature belongs to repairs the nodes
            }
            else if (support != supportType::R && value < 0)
            {
                value = 0;
                nClipped_++;
            }
            else if (support == supportType::ZeroOne && value > 1)
            {
                value = 1;
                nClipped_++;
            }

            abscissae_[nodei][scalarComponents_[si]] = value;
        }
    }

    return true;
}


Foam::scalar
Foam::multivariateMomentInversions::conditionalMeans::smallM0() const
{
    return primaryInverter_().smallM0();
}


Foam::scalar
Foam::multivariateMomentInversions::conditionalMeans::smallZeta() const
{
    return primaryInverter_().smallZeta();
}

// ************************************************************************* //
