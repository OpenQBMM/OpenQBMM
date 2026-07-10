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

#include "IlievskiEfficiency.H"
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
    defineTypeNameAndDebug(Ilievski, 0);

    addToRunTimeSelectionTable
    (
        crystalAggregationEfficiency,
        Ilievski,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Ilievski::Ilievski
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crystalAggregationEfficiency(dict, mesh),
    Utip_(readScalar(dict.lookup("Utip"))),
    psi0_(dict.lookupOrDefault<scalar>("psi0", 0.00042)),
    CIlievski_(dict.lookupOrDefault<scalar>("CIlievski", 6.25e-19)),
    rho_(dict.lookupOrDefault<scalar>("rho", 1000.0)),
    L10_(dict.lookupOrDefault<scalar>("L10", 0.0)),
    useConfiguredL10_(dict.found("L10")),
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
crystalAggregationEfficiencies::Ilievski::~Ilievski()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Ilievski::Pc
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
    const scalar G = max(GField_[celli], SMALL);
    const scalar nu = max(nu_[celli], SMALL);
    const scalar eps = max(epsilon_[celli], SMALL);

    if (!std::isfinite(G) || !std::isfinite(nu) || !std::isfinite(eps) || Utip_ <= SMALL)
    {
        return 0.0;
    }

    // Calculate shear rate from turbulent dissipation
    // gamma = sqrt(epsilon / nu)
    const scalar gamma = Foam::sqrt(eps / nu);
    if (!std::isfinite(gamma))
    {
        return 0.0;
    }

    // Dynamic viscosity mu = rho * nu
    // Assuming water density ~1000 kg/m^3
    const scalar rho = rho_;  // TODO: Get from field if available
    const scalar mu = rho * nu;
    if (!std::isfinite(mu) || mu <= SMALL)
    {
        return 0.0;
    }

    // Use local crystal size or average of d1 and d2 as L10
    const scalar L10_local =
        useConfiguredL10_
      ? max(L10_, SMALL)
      : max(0.5 * (d1 + d2), SMALL);

    if (!std::isfinite(L10_local) || L10_local <= SMALL)
    {
        return 0.0;
    }

    // Calculate Ilievski efficiency
    // Psi_10 = psi0_ / (1 + CIlievski * mu^1.5 * ((gamma^2/Utip)^2.25 / (G/L10)^3))
    
    const scalar term1 = Foam::pow(mu, 1.5);
    const scalar term2 = Foam::pow(sqr(gamma) / Utip_, 2.25);
    const scalar term3 = Foam::pow(G / L10_local, 3);

    if (!std::isfinite(term1) || !std::isfinite(term2) || !std::isfinite(term3))
    {
        return 0.0;
    }
    
    const scalar denominator = 1.0 + CIlievski_ * term1 * term2 / (term3 + SMALL);
    if (!std::isfinite(denominator) || denominator <= SMALL)
    {
        return 0.0;
    }
    
    const scalar psi = psi0_ / denominator;

    // Bound efficiency between 0 and 1
    if (!std::isfinite(psi))
    {
        return 0.0;
    }

    return max(0.0, min(1.0, psi));
}


// ************************************************************************* //
