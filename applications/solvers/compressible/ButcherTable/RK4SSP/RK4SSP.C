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

#include "RK4SSP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ButcherTables
{
    defineTypeNameAndDebug(RK4SSP, 0);
    addToRunTimeSelectionTable(ButcherTable, RK4SSP, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ButcherTables::RK4SSP::RK4SSP(const fvMesh& mesh) : ButcherTable(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ButcherTables::RK4SSP::~RK4SSP()
{}


// * * * * * * * * * * *  Protected Member Fucntions * * * * * * * * * * * * //

Foam::List<Foam::scalarList>
Foam::ButcherTables::RK4SSP::conservedVariablesCoeffs() const
{
    List<scalarList> c(10);
    c[0] = {1.0};

    for (label i = 1; i < 4; i++)
    {
        c[i] = scalarList(i + 1, 0.0);
        c[i][i] = 1.0;
    }

    c[4] = {3.0/5.0, 0.0, 0.0, 0.0, 2.0/5.0};

    for (label i = 5; i < 9; i++)
    {
        c[i] = scalarList(i + 1, 0.0);
        c[i][i] = 1.0;
    }

    c[9] = {1.0/25.0, 0, 0, 0, 9.0/25.0, 0, 0, 0, 0, 3.0/5.0};

    return c;
}

Foam::List<Foam::scalarList> Foam::ButcherTables::RK4SSP::fluxCoeffs() const
{
    List<scalarList> f(10);
    f[0] = {1.0/6.0};

    for (label i = 1; i < 4; i++)
    {
        f[i] = scalarList(i + 1, 0.0);
        f[i][i] = 1.0/6.0;
    }

    f[4] = {0.0, 0.0, 0.0, 0.0, 1.0/15.0};

    for (label i = 5; i < 9; i++)
    {
        f[i] = scalarList(i + 1, 0.0);
        f[i][i] = 1.0/6.0;
    }

    f[9] = {0, 0, 0, 0, 3.0/50.0, 0, 0 , 0, 0, 1.0/10.0};

    return f;
}

// ************************************************************************* //
