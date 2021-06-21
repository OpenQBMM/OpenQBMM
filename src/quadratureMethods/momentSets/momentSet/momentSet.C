/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
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

#include "momentSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentSet::momentSet
(
    const label nMoments,
    const label nDimensions,
    const labelListList& momentOrders,
    const word& support,
    const scalar initValue
)
:
    mappedList(nMoments, momentOrders, initValue),
    nMoments_(nMoments),
    nDimensions_(nDimensions),
    momentOrders_(momentOrders),
    support_(support)
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorInFunction
            << "The specified support is invalid." << nl
            << "    Valid supports are: R, RPlus and 01." << nl
            << "    Moment set: " << (*this)
            << abort(FatalError);
    }

    if (nDimensions_ > maxNDFDimensions_)
    {
        FatalErrorInFunction
            << "The number of maximum dimensions for the NDF is " 
            << maxNDFDimensions_ << "." << nl
            << "    Specified number of dimensions: " << nDimensions_
            << abort(FatalError);
    }
}

Foam::momentSet::momentSet
(
    const scalarList& m,
    const label nDimensions,
    const labelListList& momentOrders,
    const word& support
)
:
    mappedList(m, momentOrders),
    nMoments_(m.size()),
    nDimensions_(nDimensions),
    momentOrders_(momentOrders),
    support_(support)
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorInFunction
            << "The specified support is invalid." << nl
            << "    Valid supports are: R, RPlus and 01."
            << abort(FatalError);
    }

    if (nDimensions_ > maxNDFDimensions_)
    {
        FatalErrorInFunction
            << "The number of maximum dimensions for the NDF is " 
            << maxNDFDimensions_ << "." << nl
            << "    Specified number of dimensions: " << nDimensions_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentSet::~momentSet()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::momentSet::setSize(const label newSize)
{
    Foam::mappedList<scalar>::setSize(newSize);
    nMoments_ = newSize;
}

void Foam::momentSet::resize(const label newSize)
{
    (*this).setSize(newSize);
}

// ************************************************************************* //
