/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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

#include "environmentMixingModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::populationBalanceSubModels::environmentMixingModel>
Foam::populationBalanceSubModels::environmentMixingModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word environmentMixingModelType(dict.lookup("environmentMixingModel"));

    Info<< "Selecting environmentMixingModel "
        << environmentMixingModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(environmentMixingModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "environmentMixingModel::New(const dictionary&) : " << endl
            << "    unknown environmentMixingModelType type "
            << environmentMixingModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid environmentMixingModelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<environmentMixingModel>(cstrIter()(dict, mesh));
}


// ************************************************************************* //
