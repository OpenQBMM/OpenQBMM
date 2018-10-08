/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "noPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(noPopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        noPopulationBalance,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceModels::noPopulationBalance::noPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    populationBalanceModel(name, dict, phi),
    name_(name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceModels::noPopulationBalance::~noPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceModels::noPopulationBalance::realizableCo
() const
{
    return 1.0;
}

Foam::scalar
Foam::populationBalanceModels::noPopulationBalance::CoNum
() const
{
    return 0.0;
}

void Foam::populationBalanceModels::noPopulationBalance::solve()
{
    return;
}

// ************************************************************************* //
