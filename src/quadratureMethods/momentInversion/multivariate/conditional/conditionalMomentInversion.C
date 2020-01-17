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
    supports_(dict.lookup("supports")),
    moments_(momentOrders.size(), momentOrders, supports_[0], Zero),
    conditionalWeights_(nNodes_.size()),
    conditionalMoments_(nNodes_.size()),
    invVR_(nNodes_.size() - 1),
    momentInverters_(nNodes_.size())
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
        label mi = 0;
        Map<label> nodeMap(0);
        setNodeMap(nodeMap, nDimensions, nNodes_, 0, mi, pos);
        label nCmpts = nodeMap.size();

        conditionalWeights_.set
        (
            dimi,
            new mappedScalarList(nCmpts, nodeMap, Zero)
        );
    }

    forAll(conditionalMoments_, dimi)
    {
        nNodesCM = nNodes_;
        nNodesCM[dimi] = 2*nNodes_[dimi];

        label nDimensions = dimi + 1;
        labelList pos(nDimensions);
        label mi = 0;
        Map<label> conditionalMap(0);
        setNodeMap(conditionalMap, nDimensions, nNodesCM, 0, mi, pos);

        label nCmpts = conditionalMap.size();

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
                new mappedList<scalar>(nCmpts, conditionalMap, Zero)
            );
        }
    }

    forAll(invVR_, dimi)
    {
        label nDimensions = dimi;
        labelList pos(nDimensions);
        label mi = 0;
        Map<label> VRMap(0);
        setNodeMap(VRMap, nDimensions, nNodes_, 0, mi, pos);
        label nCmpts = VRMap.size();

        invVR_.set
        (
            dimi,
            new mappedList<scalarSquareMatrix>
            (
                nCmpts,
                VRMap,
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

    univariateMomentSet momentsToInvert(nPureMoments_[0], supports_[0]);

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
            labelList pos(dimi, 0);
            setVR(dimi - 1, pos, 0);
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

void Foam::multivariateMomentInversions::conditional::setNodeMap
(
    Map<label>& map,
    const label nDimensions,
    const labelList& nNodes,
    label dimi,
    label& mi,
    labelList& pos
)
{
    if (dimi == 0)
    {
        label size = 1;
        for (label nodei = 0; nodei < nDimensions; nodei++)
        {
            size *= nNodes[nodei];
        }

        map = Map<label>(size);
    }

    if (dimi < nDimensions)
    {
        for (label i = 0; i < nNodes[dimi]; i++)
        {
            pos[dimi] = i;
            setNodeMap(map, nDimensions, nNodes, dimi + 1, mi, pos);
        }
    }
    else
    {
        map.insert(mappedList<scalar>::listToLabel(pos), mi);
        mi++;
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

void Foam::multivariateMomentInversions::conditional::setVR
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
            setVR(dimj, pos, ai + 1);
        }
    }
    else
    {
        scalarDiagonalMatrix weights;
        scalarDiagonalMatrix x;

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

            if (mag(abscissa) > SMALL && weight > SMALL)
            {
                x.append(abscissa);
                weights.append(weight);
            }
        }

        scalarSquareMatrix invR(weights.size(), Zero);

        forAll(weights, nodei)
        {
            invR[nodei][nodei] = 1.0/weights[nodei];
        }

        Vandermonde V(x);
        scalarSquareMatrix invV(V.inv());
        labelList posVR(max(1, dimj), 0);

        if (dimj > 0)
        {
            for (label ai = 0; ai < posVR.size(); ai++)
            {
                posVR[ai] = pos[ai];
            }
        }

        invVR_[dimj](posVR) = invR*invV;
    }
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
            cycleAlphaWheeler(dimi, ai + 1, pos);
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
            supports_[dimi]
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
                weights_[nodei] = 0;
            }
        }

        return true;
    }
}


// ************************************************************************* //
