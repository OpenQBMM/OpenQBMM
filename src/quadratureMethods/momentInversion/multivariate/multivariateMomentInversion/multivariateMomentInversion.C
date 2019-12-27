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

#include "multivariateMomentInversion.H"
#include "mappedLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multivariateMomentInversion, 0);
    defineRunTimeSelectionTable(multivariateMomentInversion, dictionary);
}


// * * * * * * * * * * * * * *  Static Functions * * * * * * * * * * * * * * //

bool Foam::multivariateMomentInversion::compare
(
    const labelList& index1,
    const labelList& index2
)
{
    label size = min(index1.size(), index2.size());

    for (label i = 0; i < size; i++)
    {
        if (index1[i] != index2[i])
        {
            return false;
        }
    }

    return true;
}


bool Foam::multivariateMomentInversion::compare
(
    const labelList& index1,
    const labelList& index2,
    const label size
)
{
    for (label i = 0; i < size; i++)
    {
        if (index1[i] != index2[i])
        {
            return false;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversion::multivariateMomentInversion
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
:
    nDistributionDims_(momentOrders[0].size()),
    nvelocityDimensions_
    (
        velocityIndexes[0] == -1 ? 0 : velocityIndexes.size()
    ),
    momentOrders_(momentOrders),
    nodeIndexes_(nodeIndexes),
    velocityIndexes_(velocityIndexes),
    nNodes_(nDistributionDims_, 1),
    weights_(nodeIndexes.size(), nodeIndexes, Zero),
    abscissae_
    (
        nodeIndexes.size(),
        nodeIndexes,
        scalarList(nDistributionDims_ - nvelocityDimensions_, Zero)
    ),
    velocityAbscissae_(nodeIndexes.size(), nodeIndexes, Zero)
{
    forAll(nodeIndexes_, nodei)
    {
        forAll(nodeIndexes_[nodei], dimi)
        {
            nNodes_[dimi] = max(nNodes_[dimi], nodeIndexes[nodei][dimi] + 1);
        }
    }
    
    if (velocityIndexes_.size() == 0)
    {
        velocityIndexes_.append(-1);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversion::~multivariateMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multivariateMomentInversion::reset()
{
    forAll(weights_, nodei)
    {
        weights_[nodei] = 0.0;
        abscissae_[nodei] = scalarList(abscissae_[0].size(), Zero);
        velocityAbscissae_[nodei] = Zero;
    }
}


// ************************************************************************* //
