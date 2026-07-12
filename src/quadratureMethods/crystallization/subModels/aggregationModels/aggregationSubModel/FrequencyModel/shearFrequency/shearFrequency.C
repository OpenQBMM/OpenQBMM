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

#include "shearFrequency.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalCollisionFrequencies
{
    defineTypeNameAndDebug(shearFrequency, 0);

    addToRunTimeSelectionTable
    (
        crystalCollisionFrequency,
        shearFrequency,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
crystalCollisionFrequencies::shearFrequency::shearFrequency
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crystalCollisionFrequency(dict, mesh),
    Cshear_(dict.lookupOrDefault<scalar>("Cshear", 1.0/6.0)),
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
crystalCollisionFrequencies::shearFrequency::~shearFrequency()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::
crystalCollisionFrequencies::shearFrequency::beta
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
    const scalar nu = max(nu_[celli], SMALL);
    const scalar eps = max(epsilon_[celli], SMALL);

    if (!std::isfinite(nu) || !std::isfinite(eps))
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

    // Calculate shear-induced collision frequency
    // beta_shear = (1/6) * gamma_dot * (L1 + L2)^3
    const scalar beta = Cshear_ * gammaDot * pow3(d1 + d2);

    if (!std::isfinite(beta) || beta <= 0.0)
    {
        return 0.0;
    }

    return beta;
}


// ************************************************************************* //
