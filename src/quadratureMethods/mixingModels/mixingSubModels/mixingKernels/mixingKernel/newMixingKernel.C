/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "mixingKernel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mixingSubModels::mixingKernel>
Foam::mixingSubModels::mixingKernel::New
(
    const dictionary& dict
)
{
    word mixingKernelType(dict.lookup("mixingKernel"));

    Info<< "Selecting mixingKernel "
        << mixingKernelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(mixingKernelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "mixingKernel::New(const dictionary&) : " << endl
            << "    unknown mixingKernelType type "
            << mixingKernelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid mixingKernelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<mixingKernel>(cstrIter()(dict));
}


// ************************************************************************* //
