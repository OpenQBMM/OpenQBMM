/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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

#include "momentSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentSet::momentSet
(
    const label nMoments,
    const labelListList& momentOrders,
    const word& support,
    const scalar initValue
)
:
    scalarList(nMoments, initValue),
    nMoments_(nMoments),
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
}

Foam::momentSet::momentSet
(
    const scalarList& m,
    const labelListList& momentOrders,
    const word& support
)
:
    scalarList(m),
    nMoments_(m.size()),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentSet::~momentSet()
{}

// ************************************************************************* //
