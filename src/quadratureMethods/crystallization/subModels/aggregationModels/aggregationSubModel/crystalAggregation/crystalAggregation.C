/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Junjie Li
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

#include "crystalAggregation.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(crystalAggregation, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        crystalAggregation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::crystalAggregation::
crystalAggregation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    frequency_
    (
        crystalCollisionFrequency::New(dict, mesh)
    ),
    efficiency_
    (
        crystalAggregationEfficiency::New(dict, mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::crystalAggregation::
~crystalAggregation()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::aggregationKernels::crystalAggregation::
preUpdate()
{
    frequency_->update();
    efficiency_->update();
}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::crystalAggregation::Ka
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli,
    const label environment
) const
{
    // Ka = Ca * beta * Pc
    const scalar beta = frequency_->beta(d1, d2, Ur, celli);
    const scalar pc = efficiency_->Pc(d1, d2, Ur, celli);
    const scalar ka = Ca_.value() * beta * pc;

    if (!std::isfinite(beta) || !std::isfinite(pc) || !std::isfinite(ka) || ka <= 0.0)
    {
        return 0.0;
    }

    return ka;
}


// ************************************************************************* //
