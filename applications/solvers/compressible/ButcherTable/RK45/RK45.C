/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 Jeffrey Heylmun
    Copyright (C) 2021 Alberto Passalacqua
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

#include "RK45.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ButcherTables
{
    defineTypeNameAndDebug(RK45, 0);
    addToRunTimeSelectionTable(ButcherTable, RK45, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ButcherTables::RK45::RK45(const fvMesh& mesh) : ButcherTable(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ButcherTables::RK45::~RK45()
{}


// * * * * * * * * * * *  Protected Member Fucntions * * * * * * * * * * * * //

Foam::List<Foam::scalarList>
Foam::ButcherTables::RK45::conservedVariablesCoeffs() const
{
    return
    {
        {1.0},
        {0.0, 1.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 0.0}
    };
}

Foam::List<Foam::scalarList> Foam::ButcherTables::RK45::fluxCoeffs() const
{
    return
    {
        {0.5},
        {0.0, 0.5},
        {0.0, 0.0, 1.0},
        {1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0}
    };
}

// ************************************************************************* //
