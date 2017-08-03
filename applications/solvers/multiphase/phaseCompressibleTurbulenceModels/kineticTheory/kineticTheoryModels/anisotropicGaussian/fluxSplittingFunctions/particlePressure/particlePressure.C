/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-06-26 Jeff Heylmun:    Changed alpha to phase so that twoPhaseSystem can
                            be accessed
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

#include "particlePressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace fluxSplittingFunctions
{
    defineTypeNameAndDebug(particlePressure, 0);
    addToRunTimeSelectionTable
    (
        fluxSplittingFunction,
        particlePressure,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::fluxSplittingFunctions::particlePressure::
particlePressure
(
    const dictionary& dict
)
:
    fluxSplittingFunction(dict),
    minPpk_
    (
        dimensionedScalar::lookupOrDefault
        (
            "minPpk",
            dict,
            dimensionSet(1, -1, -2, 0, 0),
            SMALL
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::fluxSplittingFunctions::particlePressure::
~particlePressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::fluxSplittingFunctions::particlePressure::h2
(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const volScalarField& ppfr,
    const dimensionedScalar& e
) const
{
    const dimensionedScalar eta = 0.5*(1.0 + e);
    volScalarField ppk(max(alpha1*Theta*rho1, minPpk_));
    volScalarField pps(4.0*eta*alpha1*g0*ppk + ppfr);

    return pow(pps/(pps + ppk), h2Pow_);
}


// ************************************************************************* //
