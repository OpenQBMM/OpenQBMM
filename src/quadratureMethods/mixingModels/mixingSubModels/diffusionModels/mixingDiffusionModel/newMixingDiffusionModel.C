/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
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

#include "mixingDiffusionModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mixingSubModels::mixingDiffusionModel>
Foam::mixingSubModels::mixingDiffusionModel::New
(
    const dictionary& dict
)
{
    word mixingDiffusionModelType(dict.lookup("mixingDiffusionModel"));

    Info<< "Selecting mixingDiffusionModel "
        << mixingDiffusionModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(mixingDiffusionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "mixingDiffusionModel::New(const dictionary&) : " << endl
            << "    unknown mixingDiffusionModelType type "
            << mixingDiffusionModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid mixingDiffusionModelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<mixingDiffusionModel>(cstrIter()(dict));
}


// ************************************************************************* //
