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
    aggregation_(dict.lookup("aggregation")),
    breakup_(dict.lookup("breakup")),
    growth_(dict.lookup("growth")),
    nucleation_(dict.lookup("nucleation")),
    aggregationKernel_
    (
        Foam::populationBalanceSubModels::aggregationKernel::New
        (
            dict.subDict("aggregationKernel"),
            phi_.mesh()
        )
    ),
    breakupKernel_
    (
        Foam::populationBalanceSubModels::breakupKernel::New
        (
            dict.subDict("breakupKernel"),
            phi_.mesh()
        )
    ),
    growthModel_
    (
        Foam::populationBalanceSubModels::growthModel::New
        (
            dict.subDict("growthModel")
        )
    ),
    diffusionModel_
    (
        Foam::populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    nucleationModel_
    (
        Foam::populationBalanceSubModels::nucleationModel::New
        (
            dict.subDict("nucleationModel"),
            phi_.mesh()
        )
    )
{}


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
    tmp<fvScalarMatrix> implicitSource(diffusionModel_->momentDiff(moment));
    if (!collision_)
    {
        return implicitSource;
    }
    return implicitSource + collisionKernel_->implicitCollisionSource(moment);
}

void Foam::PDFTransportModels::populationBalanceModels::sizeVelocityPopulationBalance
::explicitMomentSource()
{
    if (collision_ || aggregation_ || breakup_ || growth_ || nucleation_)
    {
        odeType::solve(quadrature_, 0);
    }
    return ;
}

Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::sizeVelocityPopulationBalance::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature,
    const label environment
)
{
    scalar source = 0.0;
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
//     if (nucleation_)
//     {
//         source += nucleationModel_->nucleationSource(momentOrder[0], celli);
//     }
    if (collision_)
    {
        source += collisionKernel_->explicitCollisionSource(momentOrder, celli);
    }

    return source;
}


// ************************************************************************* //
