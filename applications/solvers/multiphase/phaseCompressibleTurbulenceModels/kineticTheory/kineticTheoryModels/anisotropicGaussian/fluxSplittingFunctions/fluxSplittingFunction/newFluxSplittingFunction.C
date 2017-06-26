/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "fluxSplittingFunction.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::kineticTheoryModels::fluxSplittingFunction>
Foam::kineticTheoryModels::fluxSplittingFunction::New
(
    const dictionary& dict
)
{
    word fluxSplittingFunctionType(dict.lookup("fluxSplittingFunction"));

    Info<< "Selecting fluxSplittingFunction "
        << fluxSplittingFunctionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fluxSplittingFunctionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "fluxSplittingFunction::New(const dictionary&) : " << endl
            << "    unknown fluxSplittingFunctionType type "
            << fluxSplittingFunctionType
            << ", constructor not in hash table" << endl << endl
            << "    Valid fluxSplittingFunctionType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<fluxSplittingFunction>(cstrIter()(dict));
}


// ************************************************************************* //
