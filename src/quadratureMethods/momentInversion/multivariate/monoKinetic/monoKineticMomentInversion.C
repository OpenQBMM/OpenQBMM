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

#include "monoKineticMomentInversion.H"
#include "mappedLists.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(monoKinetic, 0);
    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        monoKinetic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::monoKinetic::monoKinetic
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
    sizeInverter_
    (
        univariateMomentInversion::New(dict.subDict("basicQuadrature"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::monoKinetic::~monoKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::multivariateMomentInversions::monoKinetic::calcNSizeMoments
(
    const labelListList& momentOrders
)
{
    label maxOrder = 0;

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        maxOrder = max(maxOrder, momentOrder[0]);
    }

    return maxOrder + 1;
}


bool Foam::multivariateMomentInversions::monoKinetic::invert
(
    const multivariateMomentSet& moments
)
{
    reset();

    univariateMomentSet sizeMoments(nSizeMoments_, "RPlus", Zero);

    forAll(sizeMoments, mi)
    {
        sizeMoments[mi] = moments(mi);
    }

    if (!sizeMoments.isRealizable(false))
    {
        return false;
    }

    sizeInverter_->invert(sizeMoments);
    const scalarList& sizeWeights(sizeInverter_->weights());
    const scalarList& sizeAbscissae(sizeInverter_->abscissae());

    forAll(sizeWeights, nodei)
    {
        weights_[nodei] = sizeWeights[nodei];
        abscissae_[nodei][0] = sizeAbscissae[nodei];
    }

    label nSizeNodes = sizeWeights.size();

    if (nSizeNodes > 0)
    {
        scalarDiagonalMatrix x(nSizeNodes, Zero);
        scalarSquareMatrix invR(nSizeNodes, Zero);

        forAll(sizeWeights, nodei)
        {
            x[nodei] = max(sizeAbscissae[nodei], SMALL);
            invR[nodei][nodei] = 1.0/max(sizeWeights[nodei], 1e-10);
        }

        Vandermonde V(x);
        scalarSquareMatrix invVR = invR*V.inv();

        // Compute conditional velocity moments and invert
        for (label dimi = 0; dimi < nvelocityDimensions_; dimi++)
        {
            labelList pureMomentOrder(nDistributionDims_, 0);
            pureMomentOrder[dimi + 1] = 1;

            scalarRectangularMatrix M(nSizeNodes, 1, 0);
            
            for (label nodei = 0; nodei < nSizeNodes; nodei++)
            {
                pureMomentOrder[0] = nodei;
                M(nodei, 0) = moments(pureMomentOrder);
            }

            scalarRectangularMatrix nu = invVR*M;

            forAll(sizeWeights, nodei)
            {
                if (sizeWeights[nodei] > 1e-10)
                {
                    velocityAbscissae_[nodei][dimi] = nu(nodei, 0);
                }
            }
        }
    }

    return true;
}


// ************************************************************************* //
