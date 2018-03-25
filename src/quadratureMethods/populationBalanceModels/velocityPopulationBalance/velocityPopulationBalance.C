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

#include "velocityPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(velocityPopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        velocityPopulationBalance,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::velocityPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    velocityPDFTransportModel(name, dict, phi.mesh(), "R"),
    populationBalanceModel(name, dict, phi),
    name_(name),
    collision_(dict.lookup("collision")),
    collisionKernel_
    (
        Foam::populationBalanceSubModels::collisionKernel::New
        (
            dict.subDict("collisionKernel"),
            phi_.mesh(),
            quadrature_,
            dict.subDict("odeCoeffs").lookupOrDefault("solveODESource", false)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::~velocityPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::collision() const
{
    return collision_;
}

void
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::updateImplicitCollisionSource()
{
    if (!collision_)
    {
        return;
    }
    return collisionKernel_->updateFields();
}

void
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::updateExplicitCollisionSource(const label celli)
{
    if (!collision_)
    {
        return;
    }
    return collisionKernel_->updateCells(celli);
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::implicitCollisionSource
(
    const volVectorMoment& moment
)
{
    if (!collision_)
    {
        return tmp<fvScalarMatrix>
        (
            new fvScalarMatrix
            (
                moment,
                moment.dimensions()*dimVolume/dimTime
            )
        );
    }
    return collisionKernel_->implicitCollisionSource(moment);
}

Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::explicitCollisionSource
(
    const label momenti,
    const label celli
)
{
    if (!collision_)
    {
        return 0.0;
    }
    return collisionKernel_->explicitCollisionSource(momenti, celli);
}

Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::velocityPopulationBalance::cellMomentSource
(
    const label momenti,
    const label celli
)
{
    return explicitCollisionSource(momenti, celli);
}

Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::velocityPopulationBalance::realizableCo() const
{
    return velocityPDFTransportModel::realizableCo();
}

Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::velocityPopulationBalance::CoNum() const
{
    return velocityPDFTransportModel::CoNum();
}

void Foam::PDFTransportModels::populationBalanceModels
::velocityPopulationBalance::solve()
{
    velocityPDFTransportModel::solve();
}

void Foam::PDFTransportModels::populationBalanceModels
::velocityPopulationBalance::meanTransport
(
    const surfaceScalarField& phi,
    const bool wallCollisions
)
{
    velocityPDFTransportModel::meanTransport(phi, wallCollisions);
}

void Foam::PDFTransportModels::populationBalanceModels
::velocityPopulationBalance::relativeTransport
(
    const mappedPtrList<volVectorField>& Vs,
    const bool wallCollisions
)
{
    velocityPDFTransportModel::relativeTransport(Vs, wallCollisions);
}

// ************************************************************************* //
