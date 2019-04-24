/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "turbulentDiffusion.H"
#include "addToRunTimeSelectionTable.H"

#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace diffusionModels
{
    defineTypeNameAndDebug(turbulentDiffusion, 0);

    addToRunTimeSelectionTable
    (
        diffusionModel,
        turbulentDiffusion,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::diffusionModels::turbulentDiffusion
::turbulentDiffusion
(
    const dictionary& dict
)
:
    diffusionModel(dict),
    continuousPhase_(dict.lookupOrDefault("continuousPhase", word::null)),
    gammaLam_(dict.lookup("gammaLam")),
    Sc_(readScalar(dict.lookup("Sc")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::diffusionModels::turbulentDiffusion
::~turbulentDiffusion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::diffusionModels::turbulentDiffusion
::momentDiff
(
    const volScalarField& moment
) const
{
    volScalarField gamma(turbViscosity(moment)/Sc_ + gammaLam_);

    return fvm::laplacian(gamma, moment);
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::diffusionModels::turbulentDiffusion
::turbViscosity(const volScalarField& moment) const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;
    word turbName = IOobject::groupName
        (
            turbulenceModel::propertiesName,
            continuousPhase_
        );
    if (moment.mesh().foundObject<cmpTurbModel>(turbName))
    {
        const cmpTurbModel& turb =
            moment.mesh().lookupObject<cmpTurbModel>(turbName);

        return turb.mut()/turb.rho();
    }
    else if (moment.mesh().foundObject<icoTurbModel>(turbName))
    {
        const incompressible::turbulenceModel& turb =
            moment.mesh().lookupObject<icoTurbModel>(turbName);

        return turb.nut();
    }
    else
    {
        FatalErrorInFunction
            << "No valid turbulence model for turbulent diffusion calculation."
            << exit(FatalError);

        return volScalarField::null();
    }
}

// ************************************************************************* //
