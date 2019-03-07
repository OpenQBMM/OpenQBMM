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
    nNodes_(nodeIndexes.size()),
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
    weights_(nNodes_, nodeIndexes, 0.0),
    abscissae_(nNodes_, nodeIndexes, scalarList(nDistributionDims_, 0.0)),
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
{}


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

    //- Invert size moments and build VR matrix
    univariateMomentSet sizeMoments(nSizeMoments_, supports_[0], 0.0);
    labelList order(nDistributionDims_, 0);
    forAll(sizeMoments, mi)
    {
        order[0] = mi;
        sizeMoments[mi] = moments(order);
    }
    sizeInverter_->invert(sizeMoments);
    const scalarList& sizeWeights(sizeInverter_->weights());
    const scalarList& sizeAbscissae(sizeInverter_->abscissae());

    forAll(nodeIndexes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        label sizeNode = nodeIndex[0];
        if (sizeNode <= sizeInverter_->nNodes())
        {
            weights_(nodeIndex) = sizeWeights[sizeNode - 1];
            abscissae_(nodeIndex)[0] = sizeAbscissae[sizeNode - 1];
        }
    }

    scalarDiagonalMatrix x(sizeInverter_->nNodes(), 0.0);
    scalarSquareMatrix R(sizeInverter_->nNodes(), 0.0);
    scalarSquareMatrix invR(sizeInverter_->nNodes(), 0.0);

    forAll(sizeWeights, nodei)
    {
        x[nodei] = sizeAbscissae[nodei];
        R[nodei][nodei] = sizeWeights[nodei];
        invR[nodei][nodei] = 1.0/sizeWeights[nodei];
    }
    Vandermonde V(x);
    scalarSquareMatrix invVR = invR*V.inv();
    scalarSquareMatrix VR = V()*R;

    // Compute conditional velocity moments and invert
    PtrList<mappedList<scalar>> conditionalMoments(sizeInverter_->nNodes());
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

        scalarRectangularMatrix M(sizeInverter_->nNodes(), 1, 0);
        for (label sNodei = 0; sNodei < sizeInverter_->nNodes(); sNodei++)
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

            weights_(nodeIndex) *=
                velocityInverter_->weights()(velocityNodeIndex);
            for (label dimi = 0; dimi < nGeometricDimensions_; dimi++)
            {
                abscissae_(nodeIndex)[dimi + 1] =
                    velocityInverter_->abscissae()
                    (
                        velocityNodeIndex
                    ).component(dimi);
            }
        }
    }
}

void Foam::sizeCHyQMOM::reset()
{
    forAll(weights_, nodei)
    {
        weights_[nodei] = 0.0;
        abscissae_[nodei] = scalarList(nDistributionDims_, 0.0);
    }
}


// ************************************************************************* //
