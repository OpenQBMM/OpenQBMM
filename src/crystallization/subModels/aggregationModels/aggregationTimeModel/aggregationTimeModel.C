/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2023 by Your Name
    Copyright (C) 2023 Your Organization
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
#include "aggregationTimeModel.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(aggregationTimeModel, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        aggregationTimeModel,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::aggregationTimeModel::aggregationTimeModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    FField_
    (
        IOobject
        (
            "aggregationRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume/dimTime, 0.0)
    ),
    growthRate_
    (
        mesh.lookupObject<volScalarField>("growthRate")
    ),
    sigma_
    (
        mesh.lookupObject<volScalarField>("sigma")
    ),
    epsilon_
    (
        mesh.lookupObject<volScalarField>("epsilon")
    ),
    nu_
    (
        mesh.lookupObject<volScalarField>("nu")
    ),
    T_
    (
        mesh.lookupObject<volScalarField>("T")
    ),
    rhoLiquid_(dict.lookupOrDefault<scalar>("rhoLiquid", 1000.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::aggregationTimeModel::~aggregationTimeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::aggregationTimeModel::Ka
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

    const scalar eps = max(epsilon_[celli], SMALL);
    const scalar nu = max(nu_[celli], SMALL);
    const scalar T = max(T_[celli], scalar(1.0));
    const scalar mu = max(rhoLiquid_ * nu, SMALL);
    const scalar kB = 1.38e-23;

    if (!std::isfinite(eps) || !std::isfinite(nu) || !std::isfinite(T) || !std::isfinite(mu))
    {
        return 0.0;
    }

    const scalar brownian =
        (2.0 * kB * T) / (3.0 * mu)
       *sqr(1.24 * (d1 + d2))
       /(1.24 * 1.24 * d1 * d2 + SMALL);

    const scalar turbulent =
        1.91 * 0.7970 * Foam::sqrt(eps / nu) * pow3(d1 + d2);

    const scalar M = brownian + turbulent;

    if (!std::isfinite(M) || M <= 0.0)
    {
        return 0.0;
    }

    // Constants
    const scalar Mg = 1.16;
    const scalar sg = 3.08;
    const scalar betasigma_star = 0.1;

    // Calculate lambda using the larger-to-smaller ratio to keep q symmetric.
    const scalar lambda = max(d1, d2)/max(min(d1, d2), SMALL);
    if (!std::isfinite(lambda) || lambda <= SMALL)
    {
        return 0.0;
    }

    // Calculate q
    const scalar q = max(0.0, 1.0 - 0.328 * mag(Foam::log(lambda)));
    if (!std::isfinite(q))
    {
        return 0.0;
    }

    // Calculate L
    const scalar L = max(q * Foam::sqrt(d1 * d2), SMALL);

    // Calculate Gshear
    const scalar Gshear = Foam::sqrt(eps / nu);
    if (!std::isfinite(Gshear))
    {
        return 0.0;
    }

    // Reuse the local growth-rate field computed by the active growth model.
    const scalar G = max(growthRate_[celli], scalar(0.0));
    if (!std::isfinite(G))
    {
        return 0.0;
    }

    // Calculate M1
    const scalar M1 =
        (betasigma_star * eps * G)
       /(nu * sqr(2.0 * L) * sqr(Gshear) + SMALL);

    scalar alpha = 0.0;
    if (G > SMALL && sigma_[celli] > SMALL && M1 > SMALL)
    {
        alpha =
            0.5
           *(1.0 + Foam::erf((Foam::log(M1/Mg)) / (Foam::sqrt(2.0) * Foam::log(sg))));
    }

    const scalar F = Ca_.value() * M * alpha;

    if (!std::isfinite(F) || F <= 0.0)
    {
        const_cast<volScalarField&>(FField_)[celli] = 0.0;
        return 0.0;
    }
    
    // Update FField_
    const_cast<volScalarField&>(FField_)[celli] = F;

    return F;
}


// ************************************************************************* //
