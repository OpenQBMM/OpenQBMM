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
    Copyright (C) 2019-2026 Alberto Passalacqua
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

#include "conditionalMomentInversion.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(conditional, 0);
    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        conditional,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::conditional::conditional
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
    moments_
    (
        momentOrders.size(),
        momentOrders,
        supports_,
        SMALL,
        SMALL,
        Zero
    ),
    conditionalWeights_(nNodes_.size()),
    conditionalMoments_(nNodes_.size()),
    invVR_(nNodes_.size() - 1),
    momentInverters_(nNodes_.size()),
    smallM0_(SMALL),
    smallZeta_(SMALL)
{
    forAll(momentInverters_, dimi)
    {
        momentInverters_.set
        (
            dimi,
            univariateMomentInversion::New
            (
                dict.subDict("basicQuadrature" + Foam::name(dimi))
            ).ptr()
        );

        smallM0_ = max(smallM0_, momentInverters_[dimi].smallM0());
        smallZeta_ = max(smallZeta_, momentInverters_[dimi].smallZeta());

        moments_.correctSmallM0(smallM0_);
        moments_.correctSmallZeta(smallZeta_);
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

    labelList nNodesCM = nNodes_;

    forAll(nNodes_, dimi)
    {
        label nDimensions = dimi + 1;
        labelList pos(nDimensions);
        DynamicList<labelList> nodeIndexes;
        setNodeIndexes(nodeIndexes, nDimensions, nNodes_, 0, pos);
        label nCmpts = nodeIndexes.size();

        conditionalWeights_.set
        (
            dimi,
            new mappedScalarList(nCmpts, nodeIndexes, Zero)
        );
    }

    forAll(conditionalMoments_, dimi)
    {
        nNodesCM = nNodes_;
        nNodesCM[dimi] = 2*nNodes_[dimi];

        label nDimensions = dimi + 1;
        labelList pos(nDimensions);
        DynamicList<labelList> conditionalIndexes;
        setNodeIndexes(conditionalIndexes, nDimensions, nNodesCM, 0, pos);

        label nCmpts = conditionalIndexes.size();

        conditionalMoments_.set
        (
            dimi,
            new PtrList<mappedList<scalar>>(dimi + 1)
        );

        forAll(conditionalMoments_[dimi], dimj)
        {
            conditionalMoments_[dimi].set
            (
                dimj,
                new mappedList<scalar>
                (
                    nCmpts,
                    conditionalIndexes,
                    Zero
                )
            );
        }
    }

    forAll(invVR_, dimi)
    {
        label nDimensions = dimi;
        labelList pos(nDimensions);
        DynamicList<labelList> VRIndexes;
        setNodeIndexes(VRIndexes, nDimensions, nNodes_, 0, pos);
        label nCmpts = VRIndexes.size();

        invVR_.set
        (
            dimi,
            new mappedList<scalarSquareMatrix>
            (
                nCmpts,
                VRIndexes,
                scalarSquareMatrix(nNodes_[dimi], scalar(0))
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::conditional::~conditional()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multivariateMomentInversions::conditional::invert
(
    const multivariateMomentSet& moments
)
{
    reset();

    forAll(invVR_, dimi)
    {
        forAll(invVR_[dimi], ai)
        {
            invVR_[dimi][ai] = scalarSquareMatrix(nNodes_[dimi]);
        }

        forAll(conditionalMoments_[dimi], dimj)
        {
            forAll(conditionalMoments_[dimi][dimj], ai)
            {
                conditionalMoments_[dimi][dimj][ai] = 0.0;
            }
        }
    }

    forAll(conditionalWeights_, dimi)
    {
        forAll(conditionalWeights_[dimi], ai)
        {
            conditionalWeights_[dimi][ai] = 0.0;
        }
    }

    moments_ = moments;

    // Invert primary direction first
    labelList pos(nNodes_.size(), 0);

    univariateMomentSet momentsToInvert
    (
        nPureMoments_[0],
        supports_[0],
        smallM0_,
        smallZeta_
    );

    forAll(momentsToInvert, mi)
    {
        momentsToInvert[mi] = moments(mi);
    }

    if (!momentsToInvert.isRealizable(false))
    {
        return false;
    }

    momentInverters_[0].invert(momentsToInvert);

    const scalarList& weights = momentInverters_[0].weights();
    const scalarList& abscissae = momentInverters_[0].abscissae();

    vi_ = 0;
    si_ = 0;

    forAll(weights, nodei)
    {
        conditionalWeights_[0](nodei) = weights[nodei];
    }

    forAll(nodeIndexes_, nodei)
    {
        label index = nodeIndexes_[nodei][0];

        if (index < weights.size())
        {
            weights_[nodei] = weights[index];

            if (velocityIndexes_[vi_] == 0)
            {
                velocityAbscissae_[nodei][0] = abscissae[index];
            }
            else
            {
                abscissae_[nodei][0] = abscissae[index];
            }
        }
    }

    // Solve remaining directions
    for (label dimi = 1; dimi < nNodes_.size(); dimi++)
    {
        //- Set invVR matrices
        {
            labelList posVR(dimi, 0);

            if (!setVR(dimi - 1, posVR, 0))
            {
                return false;
            }
        }

        for (label dimj = 0; dimj < dimi; dimj++)
        {
            labelList posC(dimi + 1, 0);
            cycleAlphaCM(dimi, dimj, 0, posC);
        }

        if (velocityIndexes_[vi_] == dimi - 1)
        {
            vi_++;

            if (vi_ >= velocityIndexes_.size())
            {
                vi_ = 0;
            }
        }
        else
        {
            si_++;
        }

        labelList posW(dimi + 1, 0);

        if (!cycleAlphaWheeler(dimi, 0, posW))
        {
            return false;
        }
    }

    return true;
}

void Foam::multivariateMomentInversions::conditional::setNodeIndexes
(
    DynamicList<labelList>& indexes,
    const label nDimensions,
    const labelList& nNodes,
    label dimi,
    labelList& pos
)
{
    if (dimi < nDimensions)
    {
        for (label i = 0; i < nNodes[dimi]; i++)
        {
            pos[dimi] = i;
            setNodeIndexes(indexes, nDimensions, nNodes, dimi + 1, pos);
        }
    }
    else
    {
        // The node indexes used to be packed into one label, a decimal
        // digit to a direction, which a direction of ten nodes or more
        // could not be. The mapped lists take the indexes as they are.
        indexes.append(pos);
    }
}

void
Foam::multivariateMomentInversions::conditional::cycleAlphaCM
(
    const label dimi,
    const label dimj,
    label ai,
    labelList& pos
)
{
    if (dimj == ai)
    {
        cycleAlphaCM(dimi, dimj, ai+1, pos);

        return;
    }
    else if (ai < dimi)
    {
        for (label i = 0; i < nNodes_[ai]; i++)
        {
            pos[ai] = i;
            cycleAlphaCM(dimi, dimj, ai+1, pos);
        }

        return;
    }
    else if (ai == dimi)
    {
        for (label i = 0; i < nPureMoments_[dimi]; i++)
        {
            pos[dimi] = i;
            cycleAlphaCM(dimi, dimj, ai+1, pos);
        }

        return;
    }
    else
    {
        labelList posVR(max(1, dimj), 0);

        if (dimj != 0)
        {
            for (label i = 0; i < posVR.size(); i++)
            {
                posVR[i] = pos[i];
            }
        }

        const scalarSquareMatrix& invVR = invVR_[dimj](posVR);
        label size = invVR.m();
        scalarRectangularMatrix nu(size, 1, Zero);

        for (label i = 0; i < size; i++)
        {
            pos[dimj] = i;

            if (dimj == 0)
            {
                labelList posM(nNodes_.size(), 0);
                for (label mi = 0; mi < pos.size(); mi++)
                {
                    posM[mi] = pos[mi];
                }
                nu(i, 0) = moments_(posM);
            }
            else
            {
                nu(i, 0) = conditionalMoments_[dimi][dimj - 1](pos);
            }
        }

        scalarRectangularMatrix gamma = invVR*nu;

        for (label i = 0; i < nNodes_[dimj]; i++)
        {
            pos[dimj] = i;
            if (i < size)
            {
                conditionalMoments_[dimi][dimj](pos) = gamma(i, 0);
            }
            else
            {
                conditionalMoments_[dimi][dimj](pos) = 0.0;
            }
        }
    }
}

bool Foam::multivariateMomentInversions::conditional::setVR
(
    const label dimj,
    labelList& pos,
    label ai
)
{
    if (ai < dimj)
    {
        for (label i = 0; i < nNodes_[ai]; i++)
        {
            pos[ai] = i;

            if (!setVR(dimj, pos, ai + 1))
            {
                return false;
            }
        }

        return true;
    }
    else
    {
        scalarList weights(nNodes_[dimj], Zero);
        scalarList x(nNodes_[dimj], Zero);

        for (label nodei = 0; nodei < nNodes_[dimj]; nodei++)
        {
            pos[dimj] = nodei;
            scalar abscissa = 0.0;
            scalar weight = conditionalWeights_[dimj](pos);
            if (velocityIndexes_[vi_] == dimj)
            {
                abscissa = velocityAbscissae_(pos)[vi_];
            }
            else
            {
                abscissa = abscissae_(pos)[si_];
            }

            x[nodei] = abscissa;
            weights[nodei] = weight;
        }

        // The conditional solve all the quadratures share. A node carrying
        // no weight contributes nothing, and its reciprocal weight is not
        // finite, so it is left out: the weight is the only criterion, since
        // zero is an ordinary abscissa, the central node of every symmetric
        // distribution. Two nodes kept at the same abscissa are a
        // conditional quadrature that has collapsed, which the inversion
        // reports rather than solving a singular system.
        labelList solvedNodes;
        scalarSquareMatrix invVRSolved;

        if
        (
            !Vandermonde::conditionalInverse
            (
                x,
                weights,
                SMALL,
                SMALL,
                solvedNodes,
                invVRSolved
            )
        )
        {
            return false;
        }

        labelList posVR(max(1, dimj), 0);

        if (dimj > 0)
        {
            for (label posi = 0; posi < posVR.size(); posi++)
            {
                posVR[posi] = pos[posi];
            }
        }

        invVR_[dimj](posVR) = invVRSolved;
    }

    return true;
}

bool Foam::multivariateMomentInversions::conditional::cycleAlphaWheeler
(
    const label dimi,
    label ai,
    labelList& pos
)
{
    if (ai < dimi)
    {
        for (label i = 0; i < nNodes_[ai]; i++)
        {
            pos[ai] = i;

            // An unrealizable conditional moment set found deeper in the
            // recursion has to reach invert(), which reports the failed
            // inversion to the caller
            if (!cycleAlphaWheeler(dimi, ai + 1, pos))
            {
                return false;
            }
        }

        return true;
    }
    else
    {
        pos[dimi] = 0;
        scalar cm0 = conditionalMoments_[dimi][dimi - 1](pos);

        if (mag(cm0) < SMALL)
        {
            for (label nodei = 0; nodei < nNodes_[dimi]; nodei++)
            {
                pos[dimi] = nodei;
                conditionalWeights_[dimi](pos) = cm0/nNodes_[dimi];
            }

            forAll(nodeIndexes_, nodei)
            {
                const labelList& nodeIndex = nodeIndexes_[nodei];
                pos[dimi] = nodeIndex[dimi];

                if (compare(pos, nodeIndex))
                {
                    weights_[nodei] /= nNodes_[dimi];
                }
            }

            return true;
        }

        univariateMomentSet momentsToInvert
        (
            nPureMoments_[dimi],
            supports_[dimi],
            smallM0_,
            smallZeta_
        );

        forAll(momentsToInvert, mi)
        {
            pos[dimi] = mi;
            momentsToInvert[mi] = conditionalMoments_[dimi][dimi - 1](pos);
        }

        if (!momentsToInvert.isRealizable(false))
        {
            return false;
        }

        momentInverters_[dimi].invert(momentsToInvert);
        const scalarList& weights = momentInverters_[dimi].weights();
        const scalarList& abscissae = momentInverters_[dimi].abscissae();

        forAll(weights, nodei)
        {
            pos[dimi] = nodei;
            conditionalWeights_[dimi](pos) = weights[nodei];
        }

        forAll(nodeIndexes_, nodei)
        {
            const labelList& nodeIndex = nodeIndexes_[nodei];
            label index = nodeIndex[dimi];
            pos[dimi] = index;

            bool sameNode = compare(pos, nodeIndex);

            if (index < weights.size() && sameNode)
            {
                weights_[nodei] *= weights[index];

                if (velocityIndexes_[vi_] == dimi)
                {
                    velocityAbscissae_[nodei][vi_] = abscissae[index];
                }
                else
                {
                    abscissae_[nodei][si_] = abscissae[index];
                }
            }
            else if (sameNode && nodei != 0)
            {
                // The direction carries fewer quadrature nodes than the
                // case declares, so this one has nothing to be weighted
                // by. Node zero is left alone: when the direction carries
                // none at all it is the only one that keeps a weight, and
                // zeroing it too would lose what the node set holds.
                weights_[nodei] = 0;
            }
        }

        return true;
    }
}

Foam::scalar Foam::multivariateMomentInversions::conditional::smallM0() const
{
    return smallM0_;
}

Foam::scalar Foam::multivariateMomentInversions::conditional::smallZeta() const
{
    return smallZeta_;
}


// ************************************************************************* //
