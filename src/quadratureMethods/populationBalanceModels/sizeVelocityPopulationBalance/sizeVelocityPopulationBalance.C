/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 Alberto Passalacqua
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
    nucleationModel_(),
    Uc_()
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

    // The switch used to be read and the source left out, so a case that
    // asked for nucleation ran without it and did not say so
    if (nucleation_)
    {
        const dictionary& nucleationDict = dict.subDict("nucleationModel");

        nucleationModel_ =
            Foam::populationBalanceSubModels::nucleationModel::New
            (
                nucleationDict,
                phi_.mesh()
            );

        // The nuclei form out of the continuous phase and, at the size they
        // appear with, follow it without slip, so they are born with its
        // velocity: in velocity space the source is a delta at Uc. Only the
        // mean velocity is inherited, so nucleation adds no granular
        // temperature; giving the nuclei the fluctuating velocity of the
        // carrier as well would be a refinement of this closure, not
        // something it already contains.
        const word continuousPhase
        (
            nucleationDict.lookupOrDefault("continuousPhase", word::null)
        );

        const word UcName(IOobject::groupName("U", continuousPhase));
        const fvMesh& mesh = phi_.mesh();

        if (!mesh.foundObject<volVectorField>(UcName))
        {
            FatalIOErrorInFunction(nucleationDict)
                << "The nuclei are born with the velocity of the phase they "
                << "form out of, and the field " << UcName << " does not "
                << "exist." << nl
                << "    Set continuousPhase to the name of that phase." << nl
                << exit(FatalIOError);
        }

        Uc_.cref(mesh.lookupObject<volVectorField>(UcName));

        // The size the nuclei are born with comes from the nucleation model,
        // the value of any other internal coordinate is not modelled
        const volVelocityNode& node0 = quadrature_.nodes()[0];

        if (node0.sizeIndex() == -1 || node0.scalarIndexes().size() > 1)
        {
            FatalIOErrorInFunction(nucleationDict)
                << "Nucleation needs a quadrature in which the size of the "
                << "particles is the only internal coordinate other than "
                << "their velocity." << nl
                << exit(FatalIOError);
        }
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

    // Nucleation source term
    if (nucleation_)
    {
        const volVelocityNode& node0 = quadrature.nodes()[0];

        label sizeOrder = momentOrder[node0.sizeIndex()];

        // A dimensionless weight is a volume fraction, so the moments carry
        // volume and the nuclei have to be added to it in volume as well,
        // as aggregation, breakup and growth already do
        if (node0.useVolumeFraction())
        {
            sizeOrder += node0.lengthBased() ? 3 : 1;
        }

        scalar nSource =
            nucleationModel_->nucleationSource(sizeOrder, celli, environment);

        // The delta in velocity space multiplies the source of the size
        // moment by the velocity of the nuclei, raised to the order of
        // the moment in each of its components
        const vector& Uc = Uc_()[celli];
        const labelList& velocityIndexes = node0.velocityIndexes();

        forAll(velocityIndexes, cmpt)
        {
            nSource *= pow(Uc[cmpt], momentOrder[velocityIndexes[cmpt]]);
        }

        source += nSource;
    }

    return source;
}

// ************************************************************************* //
