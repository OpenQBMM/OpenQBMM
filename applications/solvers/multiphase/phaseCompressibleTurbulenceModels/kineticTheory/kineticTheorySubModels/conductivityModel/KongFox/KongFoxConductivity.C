/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "KongFoxConductivity.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
    defineTypeNameAndDebug(KongFox, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        KongFox,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::KongFox::KongFox
(
    const dictionary& dict
)
:
    conductivityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::KongFox::~KongFox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::KongFox::kappa
(
    const phaseModel& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const dimensionedScalar eta = 0.5*(1.0 + e);

    // Drag
    volScalarField beta
    (
        refCast<const twoPhaseSystem>(alpha1.fluid()).drag(alpha1).Ki()
    );
    volScalarField rTaup
    (
        "rTaup",
        max
        (
            alpha1.fluid().otherPhase(alpha1),
            alpha1.residualAlpha()
        )*beta/rho1
    );
    volScalarField rTauc
    (
        "rTauc",
        6.0*sqrt(Theta)*max(alpha1, alpha1.residualAlpha())*g0/(da*sqrtPi)
    );

    return rho1*
    (
        5.0/2.0*Theta/(3.0*rTaup + 4.0*eta*(41.0 - 33.0*eta)*rTauc)
       *(1.0 + 12.0/5.0*sqr(eta)*(4.0*eta - 3.0)*alpha1*g0)
    );

}


// ************************************************************************* //
