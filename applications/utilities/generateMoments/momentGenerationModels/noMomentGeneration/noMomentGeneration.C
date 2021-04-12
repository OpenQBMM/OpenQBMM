/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 Alberto Passalacqua
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
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

#include "noMomentGeneration.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(noMomentGeneration, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        noMomentGeneration,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::noMomentGeneration::noMomentGeneration
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    momentGenerationModel(mesh, dict, momentOrders, nNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::noMomentGeneration::~noMomentGeneration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::noMomentGeneration::updateMoments
(
    const dictionary& dict,
    const label patchi
)
{
    label size = reset(patchi);

    forAll(moments_, mi)
    {
        moments_[mi] = scalarField("moment." + Foam::name(mi), dict, size);
    }
}


void Foam::momentGenerationSubModels::noMomentGeneration::updateMoments
(
    const dictionary& dict,
    const labelList& cells
)
{
    label size = reset(cells);

    forAll(moments_, mi)
    {
        moments_[mi] = scalarField("moment." + Foam::name(mi), dict, size);
    }
}

// ************************************************************************* //
