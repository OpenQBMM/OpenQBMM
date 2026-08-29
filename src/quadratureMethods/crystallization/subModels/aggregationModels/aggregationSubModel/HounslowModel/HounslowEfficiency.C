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

#include "HounslowEfficiency.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalAggregationEfficiencies
{
    defineTypeNameAndDebug(Hounslow, 0);

    addToRunTimeSelectionTable
    (
        crystalAggregationEfficiency,
        Hounslow,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Hounslow::Hounslow
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crystalAggregationEfficiency(dict, mesh),
    Mg_(dict.lookupOrDefault<scalar>("Mg", 1.16)),
    sigmag_(dict.lookupOrDefault<scalar>("sigmag", 3.08)),
    Lstar_sigmaY_(dict.lookupOrDefault<scalar>("Lstar_sigmaY", 1.35)),
    rhoFluid_(dict.lookupOrDefault<scalar>("rhoFluid", 1000.0)),
    qBlend_(max(dict.lookupOrDefault<scalar>("qBlend", 0.05), SMALL)),
    qMin_(max(dict.lookupOrDefault<scalar>("qMin", 0.05), SMALL)),
    GField_
    (
        mesh.lookupObject<volScalarField>("growthRate")
    ),
    nu_
    (
        mesh.lookupObject<volScalarField>("nu")
    ),
    epsilon_
    (
        mesh.lookupObject<volScalarField>("epsilon")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Hounslow::~Hounslow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Hounslow::Pc
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    if (!std::isfinite(d1) || !std::isfinite(d2) || d1 <= SMALL || d2 <= SMALL)
    {
        return 0.0;
    }

    // Get local values
    const scalar G = max(GField_[celli], scalar(0.0));
    const scalar nu = max(nu_[celli], SMALL);
    const scalar eps = max(epsilon_[celli], SMALL);

    if (!std::isfinite(G) || !std::isfinite(nu) || !std::isfinite(eps) || G <= SMALL)
    {
        return 0.0;
    }

    // Calculate shear rate from turbulent dissipation
    // gamma_dot = sqrt(epsilon / nu)
    const scalar gammaDot = Foam::sqrt(eps / nu);
    if (!std::isfinite(gammaDot))
    {
        return 0.0;
    }

    // Calculate lambda (size ratio)
    const scalar lambda = max(d1, d2) / max(min(d1, d2), SMALL);
    if (!std::isfinite(lambda) || lambda <= SMALL)
    {
        return 0.0;
    }

    // Calculate q factor. Outside the Hounslow contact-length range,
    // return zero sticking. Near the boundary, apply a smooth taper and
    // evaluate L with a small q floor to avoid a near-zero-length singularity.
    const scalar q = 1.0 - 0.328 * mag(Foam::log(lambda));
    if (!std::isfinite(q) || q <= 0.0)
    {
        return 0.0;
    }

    const scalar blendFraction = min(1.0, max(0.0, q/qBlend_));
    const scalar qTaper = sqr(blendFraction)*(3.0 - 2.0*blendFraction);
    const scalar qEffective = max(q, qMin_);

    // Calculate effective contact length L
    const scalar L = qEffective * Foam::sqrt(d1 * d2);
    if (!std::isfinite(L) || L <= SMALL || !std::isfinite(qTaper))
    {
        return 0.0;
    }

    // Dynamic viscosity of the liquid phase: mu = rhoFluid * nu
    const scalar mu = rhoFluid_ * nu;
    if (!std::isfinite(mu) || mu <= SMALL)
    {
        return 0.0;
    }

    // Calculate dimensionless parameter M
    // M = (sigma_Y * G * L*) / (gamma_dot^2 * mu * L^2)
    const scalar M = (Lstar_sigmaY_ * G)
                   / (sqr(gammaDot) * mu * sqr(L) + SMALL);

    if (!std::isfinite(M) || M <= SMALL || Mg_ <= SMALL || sigmag_ <= 1.0)
    {
        return 0.0;
    }

    // Calculate efficiency using error function
    // psi = 0.5 * (1 + erf((ln(M/Mg)) / (sqrt(2) * ln(sigmag))))
    const scalar arg = Foam::log(M / Mg_) / (Foam::sqrt(2.0) * Foam::log(sigmag_));
    if (!std::isfinite(arg))
    {
        return 0.0;
    }
    const scalar psi = qTaper*0.5 * (1.0 + Foam::erf(arg));

    // Bound efficiency between 0 and 1
    if (!std::isfinite(psi))
    {
        return 0.0;
    }

    return max(0.0, min(1.0, psi));
}


// ************************************************************************* //
