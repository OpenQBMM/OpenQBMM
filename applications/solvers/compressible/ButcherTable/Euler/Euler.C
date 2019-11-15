/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 Jeff Heylmun
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

#include "Euler.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ButcherTables
{

    defineTypeNameAndDebug(Euler, 0);
    addToRunTimeSelectionTable(ButcherTable, Euler, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ButcherTables::Euler::Euler(const fvMesh& mesh)
:
    ButcherTable(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ButcherTables::Euler::~Euler()
{}


// * * * * * * * * * * *  Protected Member Fucntions * * * * * * * * * * * * //

Foam::List<Foam::scalarList>
Foam::ButcherTables::Euler::coeffs() const
{
    return {{1.0}};
}

Foam::List<Foam::scalarList>
Foam::ButcherTables::Euler::Fcoeffs() const
{
    return {{1.0}};
}

// ************************************************************************* //
