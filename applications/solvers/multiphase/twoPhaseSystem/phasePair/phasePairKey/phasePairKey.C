/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "phasePairKey.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePairKey::phasePairKey
(
    const word& name1,
    const word& name2,
    const bool ordered
)
:
    Pair<word>(name1, name2),
    ordered_(ordered)
{}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

bool Foam::operator==
(
    const phasePairKey& a,
    const phasePairKey& b
)
{
    const int cmp = Pair<word>::compare(a, b);
    return
    (
        (a.ordered() == b.ordered())
     && (a.ordered() ? (cmp == 1) : cmp)
    );
}


bool Foam::operator!=
(
    const phasePairKey& a,
    const phasePairKey& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, phasePairKey& key)
{
    const FixedList<word, 3> toks(is);

    key.first() = toks[0];
    key.second() = toks[2];
    const word& order = toks[1];

    if (order == "in")
    {
        key.ordered_ = true;
    }
    else if (order == "and")
    {
        key.ordered_ = false;
    }
    else
    {
        FatalErrorInFunction
            << "Phase pair type is not recognised. " << toks
            << "Use (phaseDispersed in phaseContinuous) for an ordered pair,"
               " or (phase1 and phase2) for an unordered pair.\n"
            << exit(FatalError);
    }

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const phasePairKey& key)
{
    os  << token::BEGIN_LIST
        << key.first()
        << token::SPACE
        << (key.ordered() ? "in" : "and")
        << token::SPACE
        << key.second()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
