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

#include "exponentialBreakup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(exponentialBreakup, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        exponentialBreakup,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::exponentialBreakup
::exponentialBreakup
(
    const dictionary& dict
)
:
    breakupKernel(dict),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", 1.0)),
    expCoeff_
    (
        dict.lookupOrDefault
        (
            "expCoeff", 
            dimensionedScalar("expC", inv(pow3(dimLength)), 0.1)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::exponentialBreakup
::~exponentialBreakup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::breakupKernels::exponentialBreakup::Kb
(
    const volScalarField& abscissa
) const
{   
    dimensionedScalar minAbs
    (
        "minAbs", 
        abscissa.dimensions(), 
        minAbscissa_.value()
    );
    
    return Cb_*pos(abscissa - minAbs)*exp(expCoeff_*pow3(abscissa));
}

// ************************************************************************* //
