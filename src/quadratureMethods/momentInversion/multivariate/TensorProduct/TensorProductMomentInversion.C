/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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
    supports_(dict.lookup("supports")),
    univariateInverters_(nNodes_.size())
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

    label vi = 0;
    label si = 0;

    forAll(univariateInverters_, dimi)
    {
        univariateMomentSet univariateMoments
        (
            nPureMoments_[dimi],
            supports_[dimi],
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

        if (max(univariateInverters_[dimi].weights()) < SMALL)
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
                    if (dimi == velocityIndexes_[vi])
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

        if (dimi == velocityIndexes_[vi])
        {
            vi++;
        }
        else
        {
            si++;
        }
    }

    if (max(nNonZeroNodes) == 0)
    {
        return true;
    }

    label totNonZeroNodes = 1;
    label nDims = 0;

    forAll(nNonZeroNodes, dimi)
    {
        if (nNonZeroNodes[dimi] > 0)
        {
            totNonZeroNodes *= nNonZeroNodes[dimi];
            nDims++;
        }
    }

    labelListList nonZeroNodeIndexes(totNonZeroNodes, labelList(nDims, 0));
    {
        label nodei = 0;
        labelList index(nDims);
        buildIndexes(nonZeroNodeIndexes, nNonZeroNodes, 0, nodei, index);
    }

    scalarList mixedMoments(nonZeroNodeIndexes.size(), Zero);
    scalarSquareMatrix R(nonZeroNodeIndexes.size(), 1.0);

    forAll(nonZeroNodeIndexes, nodei)
    {
        mixedMoments[nodei] = moments(nonZeroNodeIndexes[nodei]);
    }

    forAll(nonZeroNodeIndexes, mi)
    {
        forAll(nonZeroNodeIndexes, nodei)
        {
            vi = 0;
            si = 0;

            forAll(nonZeroNodeIndexes[nodei], dimi)
            {
                if (dimi == velocityIndexes_[vi])
                {
                    R(mi, nodei) *=
                        pow
                        (
                            velocityAbscissae_(nonZeroNodeIndexes[nodei])[vi],
                            nonZeroNodeIndexes[mi][dimi]
                        );
                    vi++;
                }
                else
                {
                    R(mi, nodei) *=
                        pow
                        (
                            abscissae_(nonZeroNodeIndexes[nodei])[si],
                            nonZeroNodeIndexes[mi][dimi]
                        );
                    si++;
                }
            }
        }
    }

    scalarList weights(nNonZeroNodes.size());
    solve(weights, R, mixedMoments);

    forAll(nonZeroNodeIndexes, nodei)
    {
        weights_(nonZeroNodeIndexes[nodei]) = weights[nodei];
    }

    return true;
}


// ************************************************************************* //
