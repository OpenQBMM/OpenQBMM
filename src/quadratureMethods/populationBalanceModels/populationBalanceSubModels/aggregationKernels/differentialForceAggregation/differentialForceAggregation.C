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

#include "differentialForceAggregation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(differentialForceAggregation, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        differentialForceAggregation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels
::differentialForceAggregation::differentialForceAggregation
(
    const dictionary& dict
)
:
    aggregationKernel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels
::differentialForceAggregation::~differentialForceAggregation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels
::aggregationKernels::differentialForceAggregation::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{
    tmp<volScalarField> aggK
        = Ca_*sqr(abscissa1 + abscissa2)*mag(sqr(abscissa1) - sqr(abscissa2));

    aggK.ref().dimensions().reset(pow3(abscissa1.dimensions())/dimTime);

    return aggK;
}

// ************************************************************************* //
