/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
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
    odeType(phi.mesh(), dict),
    collision_(dict.lookup("collision")),
    collisionKernel_
    (
        Foam::populationBalanceSubModels::collisionKernel::New
        (
            dict.subDict("collisionKernel"),
            phi_.mesh(),
            quadrature_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::~velocityPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::updateImplicitMomentSource()
{
    if (!collision_)
    {
        return;
    }

    return collisionKernel_->updateFields();
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::implicitMomentSource
(
    const volVelocityMoment& moment
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


void Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::explicitMomentSource()
{
    if (!collision_ || collisionKernel_->implicit())
    {
        return;
    }

    return odeType::solve(quadrature_, 0);
}


void
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::updateCellMomentSource(const label celli)
{
    if (!collision_)
    {
        return;
    }

    return collisionKernel_->updateCells(celli);
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::velocityPopulationBalance::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation&,
    const label
)
{
    return collisionKernel_->explicitCollisionSource(momentOrder, celli);
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


bool
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::solveMomentSources() const
{
    return odeType::solveSources_;
}


bool
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::solveMomentOde() const
{
    return odeType::solveOde_;
}


void 
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::solve()
{
    collisionKernel_->preUpdate();
    velocityPDFTransportModel::solve();
}


bool 
Foam::PDFTransportModels::populationBalanceModels::velocityPopulationBalance
::readIfModified()
{
    odeType::read
    (
        populationBalanceProperties_.subDict(type() + "Coeffs")
    );

    return true;
}


// ************************************************************************* //
