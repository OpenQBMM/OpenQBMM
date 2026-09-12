/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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

#include "constantNucleation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(constantNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        constantNucleation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::constantNucleation
::constantNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    nucleationRate_("nucleationRate", dict),
    nucleationSize_("nucleationSize", dict)
{
    // The size is written in the units of the abscissae of the case, which
    // this model cannot know, so the dimensions are read as they are given
    if (nucleationRate_.value() < 0)
    {
        FatalIOErrorInFunction(dict)
            << "The rate of nucleation cannot be negative." << nl
            << exit(FatalIOError);
    }

    if (nucleationSize_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "The particles that form have to have a size." << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::constantNucleation
::~constantNucleation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::constantNucleation
::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    return nucleationRate_.value()*pow(nucleationSize_.value(), momentOrder);
}


// ************************************************************************* //
