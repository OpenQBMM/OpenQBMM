/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 Alberto Passalacqua
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

#include "List.H"
#include "supportType.H"

namespace Foam
{

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

supportType wordToSupportType(const word& supportName)
{
    if (supportName == "R")
    {
        return supportType::R;
    }
    else if (supportName == "RPlus")
    {
        return supportType::RPlus;
    }
    else if (supportName == "ZeroOne")
    {
        return supportType::ZeroOne;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown support name: " << supportName
            << ". Valid names are: R, RPlus, ZeroOne." << endl
            << abort(FatalError);

        return supportType::RPlus;
    }
}

word supportTypeToWord(const supportType& support)
{
    switch (support)
    {
        case supportType::R:
            return "R";
        case supportType::RPlus:
            return "RPlus";
        case supportType::ZeroOne:
            return "ZeroOne";
        default:
            FatalErrorInFunction
                << "Unknown support type: " << support
                << ". Valid types are: R, RPlus, ZeroOne." << endl
                << abort(FatalError);

            return "Unknown";
    }
}

List<supportType> wordListToSupportTypeList(const wordList& supportNames)
{
    List<supportType> supports(supportNames.size());

    forAll(supportNames, i)
    {
        supports[i] = wordToSupportType(supportNames[i]);
    }

    return supports;
}

wordList supportTypeListToWordList(const List<supportType>& supports)
{
    wordList supportNames(supports.size());

    forAll(supports, i)
    {
        supportNames[i] = supportTypeToWord(supports[i]);
    }

    return supportNames;
}

Ostream& operator<<(Ostream& os, const supportType& support)
{
    os << supportTypeToWord(support);

    return os;
}

Istream& operator>>(Istream& is, supportType& support)
{
    word supportName;
    is >> supportName;
    support = wordToSupportType(supportName);

    return is;
}

} // End namespace Foam


// ************************************************************************* //
