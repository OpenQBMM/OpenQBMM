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

#include "constantAggregation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(constantAggregation, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        constantAggregation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::constantAggregation
::constantAggregation
(
    const dictionary& dict
)
:
    aggregationKernel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::constantAggregation
::~constantAggregation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::aggregationKernels::constantAggregation::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{   
    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "constantAggregationK",
                    abscissa1.mesh().time().timeName(),
                    abscissa1.mesh()
                ),
                abscissa1.mesh(),
                dimensionedScalar
                (
                    "constAggK", 
                    pow3(abscissa1.dimensions())/dimTime, 
                    Ca_.value()
                )
            )
        );
}

// ************************************************************************* //
