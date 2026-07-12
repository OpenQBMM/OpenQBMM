/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Junjie Li
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

#include "IlievskiAreaTriggeredEfficiency.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalAggregationEfficiencies
{
    defineTypeNameAndDebug(IlievskiAreaTriggered, 0);

    addToRunTimeSelectionTable
    (
        crystalAggregationEfficiency,
        IlievskiAreaTriggered,
        dictionary
    );
}
}
}
}


Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::IlievskiAreaTriggered::IlievskiAreaTriggered
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crystalAggregationEfficiency(dict, mesh),
    Utip_(readScalar(dict.lookup("Utip"))),
    L10_(dict.lookupOrDefault<scalar>("L10", 0.0)),
    useConfiguredL10_(dict.found("L10")),
    m2Name_(dict.lookupOrDefault<word>("m2Name", "moment.2.populationBalance")),
    m2Offset_(dict.lookupOrDefault<scalar>("m2Offset", 0.0)),
    alphaAreaAgg_(dict.lookupOrDefault<scalar>("alphaAreaAgg", 0.0)),
    areaPowerAgg_(dict.lookupOrDefault<scalar>("areaPowerAgg", 1.0)),
    psiBoostMax_(dict.lookupOrDefault<scalar>("psiBoostMax", 0.0)),
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
    ),
    m2Field_
    (
        mesh.lookupObject<volScalarField>(m2Name_)
    )
{}


Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::IlievskiAreaTriggered::~IlievskiAreaTriggered()
{}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::IlievskiAreaTriggered::Pc
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

    const scalar G = max(GField_[celli], SMALL);
    const scalar nu = max(nu_[celli], SMALL);
    const scalar eps = max(epsilon_[celli], SMALL);

    if (!std::isfinite(G) || !std::isfinite(nu) || !std::isfinite(eps) || Utip_ <= SMALL)
    {
        return 0.0;
    }

    const scalar gamma = Foam::sqrt(eps / nu);
    if (!std::isfinite(gamma))
    {
        return 0.0;
    }

    const scalar rho = 1000.0;
    const scalar mu = rho * nu;
    if (!std::isfinite(mu) || mu <= SMALL)
    {
        return 0.0;
    }

    const scalar L10Local =
        useConfiguredL10_
      ? max(L10_, SMALL)
      : max(0.5 * (d1 + d2), SMALL);

    if (!std::isfinite(L10Local) || L10Local <= SMALL)
    {
        return 0.0;
    }

    const scalar term1 = Foam::pow(mu, 1.5);
    const scalar term2 = Foam::pow(sqr(gamma) / Utip_, 2.25);
    const scalar term3 = Foam::pow(G / L10Local, 3);
    if (!std::isfinite(term1) || !std::isfinite(term2) || !std::isfinite(term3))
    {
        return 0.0;
    }

    const scalar denominator = 1.0 + 6.25e-19 * term1 * term2 / (term3 + SMALL);
    if (!std::isfinite(denominator) || denominator <= SMALL)
    {
        return 0.0;
    }

    const scalar psi0 = 0.00042 / denominator;
    if (!std::isfinite(psi0))
    {
        return 0.0;
    }

    const scalar m2 = max(m2Field_[celli], scalar(0.0));
    const scalar excessM2 = max(m2 - max(m2Offset_, scalar(0.0)), scalar(0.0));
    scalar activation = 0.0;
    if (alphaAreaAgg_ > SMALL && areaPowerAgg_ > SMALL && excessM2 > 0.0)
    {
        const scalar argument = Foam::pow(alphaAreaAgg_ * excessM2, areaPowerAgg_);
        if (std::isfinite(argument))
        {
            activation = 1.0 - Foam::exp(-argument);
        }
    }

    if (!std::isfinite(activation))
    {
        activation = 0.0;
    }

    const scalar psi = psi0 + max(0.0, psiBoostMax_) * max(0.0, min(1.0, activation));
    if (!std::isfinite(psi))
    {
        return 0.0;
    }

    return max(0.0, min(1.0, psi));
}


// ************************************************************************* //
