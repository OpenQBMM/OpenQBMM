/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Alberto Passalacqua
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

#include "noVelocityAdvection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{
    defineTypeNameAndDebug(noAdvection, 0);

    addToRunTimeSelectionTable
    (
        velocityMomentAdvection,
        noAdvection,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityAdvection::noAdvection::noAdvection
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const word& support
)
:
    velocityMomentAdvection(dict, quadrature, support)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityAdvection::noAdvection::~noAdvection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::velocityAdvection::noAdvection::realizableCo() const
{
    return 1.0;
}

Foam::scalar Foam::velocityAdvection::noAdvection::CoNum() const
{
    return 1.0;
}

void Foam::velocityAdvection::noAdvection::update()
{
    return;
}

void Foam::velocityAdvection::noAdvection::update
(
    const surfaceScalarField& phi,
    const bool wallCollisions
)
{
    return;
}

void Foam::velocityAdvection::noAdvection::update
(
    const mappedPtrList<volVectorField>& Us,
    const bool wallCollisions
)
{
    return;
}

// ************************************************************************* //
