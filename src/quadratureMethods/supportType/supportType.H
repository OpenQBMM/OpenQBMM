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

Class
    Foam::supportTypes

Description
    Types of support for the number density function and the associated moment
    sequences.

SourceFiles
    supportTypes.C

\*---------------------------------------------------------------------------*/

#ifndef supportType_H
#define supportType_H

#include "word.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class supportTypes Declaration
\*---------------------------------------------------------------------------*/

//- Enumeration for support types
enum class supportType
{
    R,          // Positive real line ]-∞, +∞[.
    RPlus,      // Positive real line [0, +∞[.
    ZeroOne     // Compact unit interval [0, 1].
};

//- Convert string to supportType enum
inline supportType wordToSupportType(const word& supportName);

//- Convert supportType enum to string
inline word supportTypeToWord(const supportType& support);

//- Convert wordList to supportTypeList
inline List<supportType> wordListToSupportTypeList(const wordList& supportNames);

//- Convert supportTypeList to wordList
inline wordList supportTypeListToWordList(const List<supportType>& supports);

//- Output operator for supportType
inline Ostream& operator<<(Ostream& os, const supportType& support);

//- Input operator for supportType
inline Istream& operator>>(Istream& is, supportType& support);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "supportType.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
