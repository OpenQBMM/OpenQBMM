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

#include "RKF45.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ButcherTables
{

    defineTypeNameAndDebug(RKF45, 0);
    addToRunTimeSelectionTable(ButcherTable, RKF45, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ButcherTables::RKF45::RKF45(const fvMesh& mesh)
:
    ButcherTable(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ButcherTables::RKF45::~RKF45()
{}


// * * * * * * * * * * *  Protected Member Fucntions * * * * * * * * * * * * //

Foam::List<Foam::scalarList>
Foam::ButcherTables::RKF45::conservedVariablesCoeffs() const
{
    return
    {
        {1.0},
        {0.0, 1.0},
        {0.0, 0.0, 1.0},
        {0.0, 0.0, 0.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };
}

Foam::List<Foam::scalarList> Foam::ButcherTables::RKF45::fluxCoeffs() const
{
    return
    {
        {0.25},
        {3.0/32.0, 9.0/32.0},
        {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0},
        {439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0},
        {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0},
        {16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0}
    };
}

// ************************************************************************* //
