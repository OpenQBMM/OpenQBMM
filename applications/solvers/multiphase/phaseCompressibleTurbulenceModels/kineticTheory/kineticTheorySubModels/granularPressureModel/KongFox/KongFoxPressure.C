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

#include "KongFoxPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(KongFox, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        KongFox,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::KongFox::KongFox
(
    const dictionary& dict
)
:
    granularPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::KongFox::~KongFox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::KongFox::
granularPressureCoeff
(
    const volScalarField& alpha1,
    const volScalarField& g0,
    const volScalarField& rho1,
    const dimensionedScalar& e
) const
{
    const dimensionedScalar eta = 0.5*(1.0 + e);

    if
    (
        !alpha1.mesh().foundObject<volScalarField>
        (
            IOobject::groupName("h2Fn", rho1.group())
        )
    )
    {
        FatalErrorInFunction
            << "Anisotropic Gaussian must be used with "
            << typeName_() << " model."
            << exit(FatalError);
    }

    const volScalarField& h2Fn =
        alpha1.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("h2Fn", rho1.group())
        );

    return rho1*alpha1*(h2Fn + 4.0*eta*alpha1*g0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::KongFox::
granularPressureCoeffPrime
(
    const volScalarField& alpha1,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const volScalarField& rho1,
    const dimensionedScalar& e
) const
{
    const dimensionedScalar eta = 0.5*(1.0 + e);
    if
    (
        !alpha1.mesh().foundObject<volScalarField>
        (
            IOobject::groupName("h2Fn", alpha1.group())
        )
    )
    {
        FatalErrorInFunction
            << "Anisotropic Gaussian must be used with "
            << typeName_() << " model."
            << exit(FatalError);
    }
    const volScalarField& h2Fn =
        alpha1.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("h2Fn", alpha1.group())
        );

    return rho1*(h2Fn + 4.0*alpha1*eta*(2.0*g0 + g0prime*alpha1));
}


// ************************************************************************* //
