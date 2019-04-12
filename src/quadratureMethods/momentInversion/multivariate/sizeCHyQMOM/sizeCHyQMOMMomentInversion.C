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

#include "sizeCHyQMOMMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(sizeCHyQMOM, 0);
    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        sizeCHyQMOM,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::sizeCHyQMOM::sizeCHyQMOM
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
    nSizeMoments_(calcNSizeMoments(momentOrders)),
    velocityMomentOrders_
    (
        multivariateMomentInversions::CHyQMOM::getMomentOrders
        (
            nvelocityDimensions_
        )
    ),
    nSizeNodes_(nSizeMoments_/2),
    velocityNodeIndexes_
    (
        multivariateMomentInversions::CHyQMOM::getNodeIndexes
        (
            nvelocityDimensions_
        )
    ),
    sizeInverter_
    (
        univariateMomentInversion::New(dict.subDict("basicQuadrature"))
    ),
    velocityInverter_
    (
        new multivariateMomentInversions::CHyQMOM
        (
            dict,
            velocityMomentOrders_,
            velocityNodeIndexes_,
            nvelocityDimensions_ == 1
          ? labelList({0})
          : (
                (nvelocityDimensions_ == 2)
              ? labelList({0, 1})
              : labelList({0, 1, 2})
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::sizeCHyQMOM::~sizeCHyQMOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::multivariateMomentInversions::sizeCHyQMOM::calcNSizeMoments
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


void Foam::multivariateMomentInversions::sizeCHyQMOM::invert
(
    const multivariateMomentSet& m
)
{
    reset();
    scalar m0 = m(0);
    if (m(0) < small)
    {
        forAll(weights_, nodei)
        {
            weights_[nodei] = m0/weights_.size();
        }
        return;
    }

    // Create temporary moment set and scale by m0
    multivariateMomentSet moments(m);
    forAll(moments, mi)
    {
        moments[mi] /= m0;
    }

    univariateMomentSet sizeMoments(nSizeMoments_, "RPlus", 0.0);
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
        label sizeNode = nodeIndex[0];
        if (sizeNode < sizeInverter_->nNodes())
        {
            weights_(nodeIndex) = sizeWeights[sizeNode];
            abscissae_(nodeIndex)[0] = sizeAbscissae[sizeNode];
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
            invR[nodei][nodei] = 1.0/max(sizeWeights[nodei], small);
        }
        Vandermonde V(x);
        scalarSquareMatrix invVR = invR*V.inv();

        // Compute conditional velocity moments and invert
        PtrList<mappedScalarList> conditionalMoments(nSizeNodes);
        forAll(conditionalMoments, sNodei)
        {
            conditionalMoments.set
            (
                sNodei,
                new mappedList<scalar>
                (
                    velocityMomentOrders_.size(),
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
            if (sizeWeights[sNodei] > small)
            {
                multivariateMomentSet momentsToInvert
                (
                    velocityMomentOrders_.size(),
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
                    nodeIndex[0] = sNodei;
                    for (label dimi = 1; dimi < nDistributionDims_; dimi++)
                    {
                        nodeIndex[dimi] = velocityNodeIndex[dimi - 1];
                    }

                    weights_(nodeIndex) *=
                        velocityInverter_->weights()(velocityNodeIndex);
                    velocityAbscissae_(nodeIndex) =
                        velocityInverter_->velocityAbscissae()
                        (
                            velocityNodeIndex
                        );
                }
            }
            else
            {
                forAll(velocityNodeIndexes_, nodei)
                {
                    const labelList& velocityNodeIndex = velocityNodeIndexes_[nodei];
                    labelList nodeIndex(nDistributionDims_, 0);
                    nodeIndex[0] = sNodei;
                    for (label dimi = 1; dimi < nDistributionDims_; dimi++)
                    {
                        nodeIndex[dimi] = velocityNodeIndex[dimi - 1];
                    }
                    weights_(nodeIndex) /= (weights_.size()/nSizeNodes_);
                }
            }
        }
    }
    forAll(weights_, nodei)
    {
        weights_[nodei] *= m0;
    }
}


// ************************************************************************* //
