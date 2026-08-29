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

#include "linnikovAggregation.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(linnikovAggregation, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        linnikovAggregation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::linnikovAggregation::
linnikovAggregation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    D1_(dict.lookupOrDefault<scalar>("D1", 0.1459)),
    D2_(dict.lookupOrDefault<scalar>("D2", 0.044)),
    K1D3_(dict.lookupOrDefault<scalar>("K1D3", 9.5446e-8)),
    simplified_(dict.lookupOrDefault<bool>("simplified", true)),
    sigma_
    (
        mesh.lookupObject<volScalarField>("sigma")
    )
{
    Info<< "    Linnikov aggregation model parameters:" << nl
        << "        D1 = " << D1_ << nl
        << "        D2 = " << D2_ << nl
        << "        K1D3 = " << K1D3_ << nl
        << "        simplified = " << simplified_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::linnikovAggregation::
~linnikovAggregation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::linnikovAggregation::Ka
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli,
    const label environment
) const
{
    if (!std::isfinite(d1) || !std::isfinite(d2) || d1 <= SMALL || d2 <= SMALL)
    {
        return 0.0;
    }

    // Get local supersaturation
    // sigma is supplied by the solver as the supersaturation field.
    const scalar sigma = max(sigma_[celli], SMALL);

    if (!std::isfinite(sigma))
    {
        return 0.0;
    }

    // Calculate the exponential term
    // exp(-D1 / (sigma + D2))
    const scalar expTerm = Foam::exp(-D1_ / (sigma + D2_));

    if (!std::isfinite(expTerm))
    {
        return 0.0;
    }

    scalar Ka_value;

    if (simplified_)
    {
        // Simplified first-order model:
        // Ka = K1 * D3 * exp(-D1 / (sigma + D2))
        Ka_value = K1D3_ * expTerm;
    }
    else
    {
        // Full model:
        // Ka = K1 * (1 - exp(-D3 * exp(-D1 / (sigma + D2))))
        // Note: K1D3 = K1 * D3, so we need to separate them
        // For simplicity, assume K1 = 1 and use K1D3 as D3
        Ka_value = 1.0 - Foam::exp(-K1D3_ * expTerm);
    }

    // Apply coefficient
    const scalar ka = Ca_.value() * Ka_value;

    if (!std::isfinite(ka) || ka <= 0.0)
    {
        return 0.0;
    }

    return ka;
}


// ************************************************************************* //
