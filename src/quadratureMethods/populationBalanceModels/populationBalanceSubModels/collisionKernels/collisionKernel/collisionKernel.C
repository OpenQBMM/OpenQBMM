/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
     \\/     M anipulation  |
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

#include "collisionKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(collisionKernel, 0);

    defineRunTimeSelectionTable(collisionKernel, dictionary);
}
}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::populationBalanceSubModels::collisionKernel::meanVelocity
(
    const scalar& m0,
    const scalar& m1
) const
{
    return m1/m0;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::collisionKernel::meanVelocity
(
    const volScalarField& m0,
    const volScalarField& m1
)
{
    return m1/m0;
}

Foam::scalar Foam::populationBalanceSubModels::collisionKernel::variance
(
    const scalar& m0,
    const scalar& sqrM1,
    const scalar& m2
) const
{
    return max(m2/m0 - sqrM1, 0.0);
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::collisionKernel::variance
(
    const volScalarField& m0,
    const volScalarField& sqrM1,
    const volScalarField& m2
)
{
    return max(m2/m0 - sqrM1, dimensionedScalar("zero", sqr(dimVelocity), 0.0));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernel::collisionKernel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature,
    const bool ode
)
:
    dict_(dict),
    mesh_(mesh),
    quadrature_(quadrature),
    momentOrders_(quadrature.momentOrders()),
    nDimensions_(momentOrders_[0].size())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernel::~collisionKernel()
{}


// ************************************************************************* //
