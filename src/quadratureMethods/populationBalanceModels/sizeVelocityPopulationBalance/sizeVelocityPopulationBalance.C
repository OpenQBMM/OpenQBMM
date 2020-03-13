/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 Alberto Passalacqua
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

#include "sizeVelocityPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(sizeVelocityPopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        sizeVelocityPopulationBalance,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::sizeVelocityPopulationBalance
::sizeVelocityPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    velocityPopulationBalance(name, dict, phi),
    aggregation_(dict.lookupOrDefault("aggregation", false)),
    breakup_(dict.lookupOrDefault("breakup", false)),
    growth_(dict.lookupOrDefault("growth", false)),
    nucleation_(dict.lookupOrDefault("nucleation", false)),
    aggregationKernel_(),
    breakupKernel_(),
    growthModel_(),
    nucleationModel_()
{
    if (aggregation_)
    {
        aggregationKernel_ =
            Foam::populationBalanceSubModels::aggregationKernel::New
            (
                dict.subDict("aggregationKernel"),
                phi_.mesh()
            );
    }

    if (breakup_)
    {
        breakupKernel_ =
            Foam::populationBalanceSubModels::breakupKernel::New
            (
                dict.subDict("breakupKernel"),
                phi_.mesh()
            );
    }

    if (growth_)
    {
        growthModel_ =
            Foam::populationBalanceSubModels::growthModel::New
            (
                dict.subDict("growthModel"),
                phi_.mesh()
            );
    }

    if (dict.found("diffusionModel"))
    {
        diffusionModel_ =
            Foam::populationBalanceSubModels::diffusionModel::New
            (
                dict.subDict("diffusionModel")
            );
    }

    if (nucleation_)
    {
        nucleationModel_ =
            Foam::populationBalanceSubModels::nucleationModel::New
            (
                dict.subDict("nucleationModel"),
                phi_.mesh()
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::sizeVelocityPopulationBalance
::~sizeVelocityPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::sizeVelocityPopulationBalance
::implicitMomentSource
(
    const volVelocityMoment& moment
)
{
    tmp<fvScalarMatrix> momentEqn
    (
        velocityPopulationBalance::implicitMomentSource(moment)
    );

    if (diffusionModel_.valid())
    {
        return momentEqn + diffusionModel_->momentDiff(moment);
    }
    else 
    {
        return momentEqn;
    }
}

void 
Foam::PDFTransportModels::populationBalanceModels::sizeVelocityPopulationBalance
::explicitMomentSource()
{
    if
    (
        (collision_ && !collisionKernel_->implicit())
      || aggregation_ || breakup_ || growth_ || nucleation_
    )
    {
        odeType::solve(quadrature_, 0);
    }

    return;
}

Foam::scalar 
Foam::PDFTransportModels::populationBalanceModels
::sizeVelocityPopulationBalance::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature,
    const label environment
)
{
    scalar source(0);

//     if (nucleation_)
//     {
//         source += nucleationModel_->nucleationSource(momentOrder[0], celli);
//     }

    // Collision source term
    if (collision_)
    {
        source += collisionKernel_->explicitCollisionSource(momentOrder, celli);
    }

    // Aggregation source term
    if (aggregation_)
    {
        source +=
            aggregationKernel_->aggregationSource
            (
                momentOrder,
                celli,
                quadrature,
                environment
            );
    }

    // Breaku source term
    if (breakup_)
    {
        source +=
            breakupKernel_->breakupSource
            (
                momentOrder,
                celli,
                quadrature
            );
    }

    // Phase space convection/growth source term
    if (growth_)
    {
        source +=
            growthModel_->phaseSpaceConvection
            (
                momentOrder,
                celli,
                quadrature
            );
    }

    return source;
}

// ************************************************************************* //
