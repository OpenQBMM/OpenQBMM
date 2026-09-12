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

#include "sizeCHyQMOMMomentInversion.H"
#include "supportType.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class velocityInversion>
Foam::multivariateMomentInversions::sizeCHyQMOMBase<velocityInversion>::
sizeCHyQMOMBase
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
        velocityInversion::getMomentOrders
        (
            nVelocityDimensions_
        )
    ),
    nSizeNodes_(nSizeMoments_/2),
    velocityNodeIndexes_
    (
        velocityInversion::getNodeIndexes
        (
            nVelocityDimensions_
        )
    ),
    sizeInverter_
    (
        univariateMomentInversion::New(sizeQuadratureDict(dict))
    ),
    velocityInverter_
    (
        new velocityInversion
        (
            dict,
            velocityMomentOrders_,
            velocityNodeIndexes_,
            nVelocityDimensions_ == 1
          ? labelList({0})
          : (
                (nVelocityDimensions_ == 2)
              ? labelList({0, 1})
              : labelList({0, 1, 2})
            )
        )
    ),
    // Both thresholds are those of the quadrature of the size direction,
    // which is the one the zero-order moment and the zeta_k of the whole
    // distribution belong to, and sizeQuadratureDict has already given it
    // the value written where the two quadratures meet.
    //
    // They used to be the larger of the values of the two sub-inverters,
    // which meant a value written in either place could only ever raise
    // the threshold: the default of the other one silently overrode it.
    smallM0_(sizeInverter_().smallM0()),
    // The zeta_k carry the dimensions of the abscissa, so a floor on them
    // is a threshold on a size. It is a threshold on a size written in
    // units of its own mean, which invert normalises the size moments to,
    // and so means the same thing whatever the size of a particle is
    // written in.
    smallZeta_(max(SMALL, sizeInverter_().smallZeta()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class velocityInversion>
Foam::multivariateMomentInversions::sizeCHyQMOMBase<velocityInversion>::
~sizeCHyQMOMBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class velocityInversion>
Foam::label Foam::multivariateMomentInversions::sizeCHyQMOMBase<velocityInversion>::
calcNSizeMoments
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

    // A univariate quadrature of n nodes is built from 2n moments, so an
    // odd count leaves a moment out of the closure. It is refused rather
    // than dropped by the integer division that sets nSizeNodes_.
    if ((maxOrder + 1) % 2 != 0)
    {
        FatalErrorInFunction
            << "The size direction carries " << maxOrder + 1
            << " moments, of orders zero to " << maxOrder << "." << nl
            << "    An even number is needed: a quadrature of n nodes is "
            << "built from 2n moments." << nl
            << exit(FatalError);
    }

    return maxOrder + 1;
}


template<class velocityInversion>
Foam::dictionary
Foam::multivariateMomentInversions::sizeCHyQMOMBase<velocityInversion>::
sizeQuadratureDict
(
    const dictionary& dict
)
{
    dictionary sizeDict(dict.subDict("basicQuadrature"));

    // A threshold written where the size and the velocity quadratures meet
    // governs the inversion of the size direction as well, so that one
    // entry sets it for the whole of the inversion: the quadrature of the
    // size direction keeps a value of its own, and cuts the zero-order
    // moment at it before anything above it is reached. Written in neither
    // place, each keeps the default it was given.
    for (const word& key : {word("smallM0"), word("smallZeta")})
    {
        if (dict.found(key))
        {
            sizeDict.set(key, dict.get<scalar>(key));
        }
    }

    return sizeDict;
}


template<class velocityInversion>
bool Foam::multivariateMomentInversions::sizeCHyQMOMBase<velocityInversion>::
invert
(
    const multivariateMomentSet& moments
)
{
    reset();
    scalar m0 = moments(0);

    // A negligible zero-order moment is dropped rather than carried by
    // nodes of null size and velocity: in a volume-fraction formulation the
    // number density of such a node is its weight over the cube of a null
    // diameter, so the mass parked there is not negligible for the
    // collision and aggregation sources of the cell it is advected into.
    if (mag(m0) < smallM0())
    {
        return true;
    }

    // The size direction is inverted in units of its own mean, and the
    // abscissae are written back in the units of the case below.
    //
    // Whether a moment set is realizable is a property of the distribution
    // and not of the units its abscissa is measured in, but the zeta_k the
    // check compares with smallZeta carry those units: a coordinate small
    // enough in them, as a volume of a micron-sized particle is in cubic
    // metres, is declared degenerate however well spread it is. Dividing
    // the moment of order k by the mean size to the power k leaves zeta_0
    // at one and the rest of order one, so that the threshold means the
    // same thing whether the size is a length, a volume or a mass.
    //
    // A mean size that is not positive is left alone: the moment set is
    // then unrealizable over R+, which the check below is what says.
    const scalar meanSize = moments(1)/m0;
    const scalar sizeScale = meanSize > 0 ? meanSize : 1.0;

    univariateMomentSet sizeMoments
    (
        nSizeMoments_,
        supportType::RPlus,
        smallM0(),
        smallZeta(),
        Zero);

    scalar scaleToPower = 1.0;

    forAll(sizeMoments, mi)
    {
        sizeMoments[mi] = moments(mi)/scaleToPower;
        scaleToPower *= sizeScale;
    }

    if (!sizeMoments.isRealizable(false))
    {
        return false;
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
            abscissae_(nodeIndex)[0] = sizeAbscissae[sizeNode]*sizeScale;
        }
    }

    label nSizeNodes = sizeWeights.size();

    if (nSizeNodes > 0)
    {
        scalarDiagonalMatrix x(nSizeNodes, Zero);
        scalarSquareMatrix invR(nSizeNodes, Zero);

        // The system is built and solved in the units the size direction
        // was inverted in, where the abscissae are of order one whatever
        // the size of a particle is written in. The conditional velocity
        // moments it returns are the same either way, since the moments it
        // is given are divided by the same powers of the scale.
        forAll(sizeWeights, nodei)
        {
            x[nodei] = sizeAbscissae[nodei];

            // The row of a size node that carries no particles is not read
            // back below, so it only has to be finite. Everywhere else the
            // weight is divided by as it is: flooring it scaled the
            // conditional velocity moments of a node that does carry
            // particles by whatever the ratio of the two happened to be,
            // which a weight is free to be below wherever the measure is
            // small in its own units.
            invR[nodei][nodei] =
                sizeWeights[nodei] > smallM0()
              ? 1.0/sizeWeights[nodei]
              : 0.0;
        }

        // The Vandermonde system has to be built from the abscissae the
        // quadrature actually carries, which are the ones stored above.
        // Moving them away from zero would both bias the conditional
        // velocity moments recovered from it and, when two of them are
        // small, make the system singular by mapping them onto the same
        // value.
        //
        // The abscissae a successful univariate inversion returns are
        // distinct, so this only refuses a size quadrature that has
        // collapsed without the realizability check having reduced it.
        if (!distinctAbscissae(x))
        {
            return false;
        }

        Vandermonde V(x);
        scalarSquareMatrix invVR = invR*V.invert();

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
                    Zero
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

            scalar scaleToPowerOfNode = 1.0;

            for (label sNodei = 0; sNodei < nSizeNodes; sNodei++)
            {
                pureMomentOrder[0] = sNodei;
                M(sNodei, 0) = moments(pureMomentOrder)/scaleToPowerOfNode;
                scaleToPowerOfNode *= sizeScale;
            }

            scalarRectangularMatrix nu = invVR*M;

            forAll(conditionalMoments, sNodei)
            {
                conditionalMoments[sNodei](velocityMomentOrder) = nu(sNodei, 0);
            }
        }

        forAll(conditionalMoments, sNodei)
        {
            if (sizeWeights[sNodei] > smallM0())
            {
                multivariateMomentSet momentsToInvert
                (
                    velocityMomentOrders_.size(),
                    velocityMomentOrders_,
                    List<supportType>(nVelocityDimensions_, supportType::R),
                    smallM0(),
                    smallZeta()
                );

                forAll(momentsToInvert, mi)
                {
                    momentsToInvert(velocityMomentOrders_[mi]) =
                        conditionalMoments[sNodei](velocityMomentOrders_[mi]);
                }

                velocityInverter_->invert(momentsToInvert);

                forAll(velocityNodeIndexes_, nodei)
                {
                    const labelList& velocityNodeIndex =
                        velocityNodeIndexes_[nodei];

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
                // A size node of negligible weight is dropped for the same
                // reason: its velocity moments cannot be inverted, and
                // keeping its weight on nodes of null velocity would park
                // mass at rest.
                forAll(velocityNodeIndexes_, nodei)
                {
                    const labelList& velocityNodeIndex =
                        velocityNodeIndexes_[nodei];

                    labelList nodeIndex(nDistributionDims_, 0);
                    nodeIndex[0] = sNodei;

                    for (label dimi = 1; dimi < nDistributionDims_; dimi++)
                    {
                        nodeIndex[dimi] = velocityNodeIndex[dimi - 1];
                    }

                    weights_(nodeIndex) = Zero;
                }
            }
        }
    }
    else
    {
        forAll(weights_, nodei)
        {
            weights_[nodei] = m0/weights_.size();
        }
    }

    return true;
}

template<class velocityInversion>
Foam::scalar
Foam::multivariateMomentInversions::sizeCHyQMOMBase<velocityInversion>::smallM0()
const
{
    return smallM0_;
}

template<class velocityInversion>
Foam::scalar
Foam::multivariateMomentInversions::sizeCHyQMOMBase<velocityInversion>::smallZeta()
const
{
    return smallZeta_;
}

// ************************************************************************* //
