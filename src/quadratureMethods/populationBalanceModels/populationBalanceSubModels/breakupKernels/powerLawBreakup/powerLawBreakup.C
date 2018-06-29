/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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

#include "powerLawBreakup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(powerLawBreakup, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        powerLawBreakup,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::powerLawBreakup
::powerLawBreakup
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    breakupKernel(dict, mesh),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", 1.0)),
    abscissaExponent_(dict.lookupOrDefault("abscissaExponent", 3.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::powerLawBreakup
::~powerLawBreakup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::breakupKernels::powerLawBreakup::Kb
(
    const scalar& abscissa,
    const label celli,
    const label environment
) const
{
    return Cb_.value()*pos0(abscissa - minAbscissa_.value())
            *pow(abscissa, abscissaExponent_);
}

// ************************************************************************* //
