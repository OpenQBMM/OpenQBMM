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

#include "fieldMomentInversion.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fieldMomentInversion> Foam::fieldMomentInversion::New
(
    const dictionary& dict,
    const label nMoments,
    const label nSecondaryNodes
)
{
    word fieldMomentInversionType
    (
        dict.lookup("fieldMomentInversion")
    );

    Info<< "Selecting fieldMomentInversion: "
        << fieldMomentInversionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fieldMomentInversionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fieldMomentInversion type "
            << fieldMomentInversionType << endl << endl
            << "Valid fieldMomentInversion types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, nMoments, nSecondaryNodes);
}


// ************************************************************************* //
