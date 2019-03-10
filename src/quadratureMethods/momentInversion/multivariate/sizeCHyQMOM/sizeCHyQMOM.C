/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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

#include "sizeCHyQMOM.H"
#include "mappedLists.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sizeCHyQMOM::sizeCHyQMOM
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes
)
:
    nDistributionDims_(momentOrders[0].size()),
    nGeometricDimensions_(nDistributionDims_ - 1),
    nMoments_(momentOrders.size()),
    nSizeMoments_(calcNSizeMoments(momentOrders)),
    nVelocityMoments_
    (
        hyperbolicConditionalMomentInversion::getNMoments
        (
            nGeometricDimensions_
        )
    ),
    momentOrders_(momentOrders),
    velocityMomentOrders_
    (
        hyperbolicConditionalMomentInversion::getMomentOrders
        (
            nGeometricDimensions_
        )
    ),
    nNodes_(nodeIndexes[0].size()),
    nSizeNodes_(nSizeMoments_/2),
    nodeIndexes_(nodeIndexes),
    velocityNodeIndexes_
    (
        hyperbolicConditionalMomentInversion::getNodeIndexes
        (
            nGeometricDimensions_
        )
    ),
    supports_({"RPlus", "R", "R", "R"}),
    weights_(nodeIndexes.size(), nodeIndexes, 0.0),
    sizeAbscissae_(nodeIndexes.size(), nodeIndexes, 0.0),
    velocityAbscissae_(nodeIndexes.size(), nodeIndexes, Zero),
    sizeInverter_
    (
        univariateMomentInversion::New(dict.subDict("basicQuadrature"))
    ),
    velocityInverter_
    (
        new hyperbolicConditionalMomentInversion
        (
            dict,
            nGeometricDimensions_
        )
    )
{
    forAll(nodeIndexes_, nodei)
    {
        forAll(nNodes_, dimi)
        {
            nNodes_[dimi] = max(nNodes_[dimi], nodeIndexes[nodei][dimi]);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sizeCHyQMOM::~sizeCHyQMOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::sizeCHyQMOM::calcNSizeMoments
(
    const labelListList& momentOrders
)
{
    label maxOrder = 0;
    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        if (momentOrder[0] > maxOrder)
        {
            maxOrder = momentOrder[0];
        }
    }
    return maxOrder + 1;
}


void Foam::sizeCHyQMOM::invert
(
    const multivariateMomentSet& moments
)
{
    reset();

    univariateMomentSet sizeMoments(nSizeMoments_, supports_[0], 0.0);
    forAll(sizeMoments, mi)
    {
        sizeMoments[mi] = moments(mi);
    }

    sizeInverter_->invert(sizeMoments);
    const scalarList& sizeWeights(sizeInverter_->weights());
    const scalarList& sizeAbscissae(sizeInverter_->abscissae());

    forAll(nodeIndexes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        label sizeNode = nodeIndex[0] - 1;
        if (sizeNode < sizeInverter_->nNodes())
        {
            weights_(nodeIndex) = sizeWeights[sizeNode];
            sizeAbscissae_(nodeIndex) = sizeAbscissae[sizeNode];
        }
    }
    label nSizeNodes = sizeWeights.size();

    if (nSizeNodes > 0)
    {
        scalarDiagonalMatrix x(nSizeNodes, 0.0);
        scalarSquareMatrix invR(nSizeNodes, 0.0);

        forAll(sizeWeights, nodei)
        {
            x[nodei] = max(sizeAbscissae[nodei], small);
            invR[nodei][nodei] = 1.0/max(sizeWeights[nodei], 1e-10);
        }
        Vandermonde V(x);
        scalarSquareMatrix invVR = invR*V.inv();

        // Compute conditional velocity moments and invert
        PtrList<mappedList<scalar>> conditionalMoments(nSizeNodes);
        forAll(conditionalMoments, sNodei)
        {
            conditionalMoments.set
            (
                sNodei,
                new mappedList<scalar>
                (
                    nVelocityMoments_,
                    velocityMomentOrders_,
                    0.0
                )
            );
        }

        forAll(velocityMomentOrders_, mi)
        {
            const labelList& velocityMomentOrder = velocityMomentOrders_[mi];
            labelList pureMomentOrder(nDistributionDims_, 0);
            for (label dimi = 1; dimi < nDistributionDims_; dimi++)
            {
                pureMomentOrder[dimi] = velocityMomentOrder[dimi - 1];
            }

            scalarRectangularMatrix M(nSizeNodes, 1, 0);
            for (label sNodei = 0; sNodei < nSizeNodes; sNodei++)
            {
                pureMomentOrder[0] = sNodei;
                M(sNodei, 0) = moments(pureMomentOrder);
            }
            scalarRectangularMatrix nu = invVR*M;

            forAll(conditionalMoments, sNodei)
            {
                conditionalMoments[sNodei](velocityMomentOrder) = nu(sNodei, 0);
            }
        }

        forAll(conditionalMoments, sNodei)
        {
            scalar conditionalM0 =  conditionalMoments[sNodei](0);
            forAll(velocityMomentOrders_, mi)
            {
                conditionalMoments[sNodei][mi] /= max(conditionalM0, 1e-10);
            }
            conditionalMoments[sNodei](0) = 1.0;
        }

        forAll(conditionalMoments, sNodei)
        {
            multivariateMomentSet momentsToInvert
            (
                nVelocityMoments_,
                velocityMomentOrders_,
                "R"
            );
            forAll(momentsToInvert, mi)
            {
                momentsToInvert(velocityMomentOrders_[mi]) =
                    conditionalMoments[sNodei](velocityMomentOrders_[mi]);
            }
            velocityInverter_->invert(momentsToInvert);

            forAll(velocityNodeIndexes_, nodei)
            {
                const labelList& velocityNodeIndex = velocityNodeIndexes_[nodei];
                labelList nodeIndex(nDistributionDims_, 0);
                nodeIndex[0] = sNodei + 1;
                for (label dimi = 1; dimi < nDistributionDims_; dimi++)
                {
                    nodeIndex[dimi] = velocityNodeIndex[dimi - 1];
                }

                if (sizeWeights[sNodei] < 1e-10)
                {
                    weights_(nodeIndex) /= (weights_.size()/nSizeNodes_);
                    velocityAbscissae_(nodeIndex) = Zero;
                }
                else
                {
                    weights_(nodeIndex) *=
                        velocityInverter_->weights()(velocityNodeIndex);
                    velocityAbscissae_(nodeIndex) =
                        velocityInverter_->abscissae()(velocityNodeIndex);
                }
            }
        }
    }
}

void Foam::sizeCHyQMOM::reset()
{
    forAll(weights_, nodei)
    {
        weights_[nodei] = 0.0;
        sizeAbscissae_[nodei] = 0.0;
        velocityAbscissae_[nodei] = Zero;
    }
}


// ************************************************************************* //
