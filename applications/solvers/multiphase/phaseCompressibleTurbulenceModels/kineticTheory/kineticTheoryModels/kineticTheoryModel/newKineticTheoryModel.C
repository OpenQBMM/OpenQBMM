/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "kineticTheoryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::kineticTheoryModel>
Foam::kineticTheoryModel::New
(
    const dictionary& dict,
    const phaseModel& phase
)
{
    word kineticTheoryModelType(dict.lookup("kineticTheoryModel"));

    Info<< "Selecting kineticTheoryModel "
        << kineticTheoryModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(kineticTheoryModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "kineticTheoryModel::New(const dictionary&) : " << endl
            << "    unknown kineticTheoryModelType type "
            << kineticTheoryModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid kineticTheoryModelType types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return
        autoPtr<kineticTheoryModel>(cstrIter()(dict, phase));
}


// ************************************************************************* //
