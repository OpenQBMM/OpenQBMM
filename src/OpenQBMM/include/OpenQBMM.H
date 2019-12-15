/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Alberto Passalacqua
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

Description
    General macros and definitions, especially to help when dealing
    with different OpenFOAM variants and versions.

SourceFiles
    OpenQBMM.H

\*---------------------------------------------------------------------------*/

#ifndef OpenQBMM_H
#define OpenQBMM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       For openfoam.com (1712 and newer)
\*---------------------------------------------------------------------------*/

#if OPENFOAM

#if (OPENFOAM >= 1806)
#define DICT_GET(Dict,Type,Key)            Dict.get<Type>(Key)
#endif

#if (OPENFOAM >= 1812)
// Construct dimensioned type with given name from dictionary lookup
#define DICT_DIMTYPE(Dict,Key)              Key,Dict
#endif

#if (OPENFOAM > 1806)
#define DICT_READ_ENTRY(Dict,Key,Val)       Dict.readEntry(Key, Val)
#endif

#if (OPENFOAM > 1606)
#define WRITE_ENTRY(Os,Key,Val)             Os.writeEntry(Key, Val)
#endif


// Housekeeping

#if (OPENFOAM < 1812)
#define DICT_DIMTYPE(Dict,Key)              Dict.lookup(Key)
#endif

#if (OPENFOAM == 1806)
// This was a bad naming choice (compiler complaints).
// Subsequently changed to readEntry(). Mnemonic: reverse of writeEntry()
#define DICT_READ_ENTRY(Dict,Key,Val)       Dict.read(Key, Val)
#elif (OPENFOAM < 1806)
#define DICT_GET(Dict,Type,Key)             Dict.lookupType<Type>(Key)
#define DICT_READ_ENTRY(Dict,Key,Val)       Dict.lookup(Key)() >> Val
#endif


/*---------------------------------------------------------------------------*\
                               For openfoam.org
\*---------------------------------------------------------------------------*/

#else

// Until 2019-11-27
#define DICT_GET(Dict,Type,Key)             Dict.lookupType<Type>(Key)

// Since 2019-11-27
///#define DICT_GET(Dict,Type,Key)             Dict.lookup<Type>(Key)

// Construct dimensioned type from dictionary lookup
#define DICT_DIMTYPE(Dict, Key)              Dict.lookup(Key)

// Compat with openfoam.com
#define DICT_READ_ENTRY(Dict,Key,Val)       Dict.lookup(Key)() >> Val

// Until 2019-07-09
#define WRITE_ENTRY(os,Key,Val)   os.writeKeyword(Key) << (Val) << ';' << nl

// Since 2019-07-09
///#define WRITE_ENTRY(os,Key,Val)             writeEntry(os, Key, Val)

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
