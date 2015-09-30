/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Matteo Icardi
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

#include "constantGrowth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(constantGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        constantGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::constantGrowth
::constantGrowth
(
    const dictionary& dict
)
:
    growthModel(dict),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", 1e-10)),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", 1e10))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::constantGrowth
::~constantGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::constantGrowth::Kg
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

    dimensionedScalar maxAbs
    (
        "maxAbs", 
        abscissa.dimensions(), 
        maxAbscissa_.value()
    );

    dimensionedScalar CgND
    (
        "CgND", 
        inv(abscissa.dimensions())*inv(abscissa.dimensions())*Cg_.dimensions(), 
        Cg_.value()
    );

    return CgND*pos(-abscissa+maxAbs)*pos(abscissa-minAbs);
}

// ************************************************************************* //
