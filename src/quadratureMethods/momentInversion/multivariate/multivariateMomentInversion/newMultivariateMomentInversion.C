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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multivariateMomentInversion>
Foam::multivariateMomentInversion::New
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
{
    word multivariateMomentInversionType
    (
        dict.lookup("type")
    );

    Info<< "Selecting multivariateMomentInversion "
        << multivariateMomentInversionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(multivariateMomentInversionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown multivariateMomentInversionType type "
            << multivariateMomentInversionType << endl << endl
            << "Valid multivariateMomentInversionType types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return
        autoPtr<multivariateMomentInversion>
        (
            cstrIter()
            (
                dict,
                momentOrders,
                nodeIndexes,
                velocityIndexes
            )
        );
}


// ************************************************************************* //
