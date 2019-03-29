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

#include "multivariateMomentInversion.H"
#include "mappedLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multivariateMomentInversion, 0);
    defineRunTimeSelectionTable(multivariateMomentInversion, dictionary);
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
    nGeometricDimensions_(velocityIndexes.size()),
    momentOrders_(momentOrders),
    nodeIndexes_(nodeIndexes),
    velocityIndexes_(velocityIndexes),
    nNodes_(nDistributionDims_, 1),
    weights_(nodeIndexes.size(), nodeIndexes, 0.0),
    abscissae_
    (
        nodeIndexes.size(),
        nodeIndexes,
        scalarList(nDistributionDims_ - nGeometricDimensions_, 0.0)
    ),
    velocityAbscissae_(nodeIndexes.size(), nodeIndexes, Zero)
{
    forAll(nodeIndexes_, nodei)
    {
        forAll(nNodes_, dimi)
        {
            nNodes_[dimi] = max(nNodes_[dimi], nodeIndexes[nodei][dimi]);
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
        abscissae_[nodei] = scalarList(abscissae_[0].size(), 0.0);
        velocityAbscissae_[nodei] = Zero;
    }
}


// ************************************************************************* //
