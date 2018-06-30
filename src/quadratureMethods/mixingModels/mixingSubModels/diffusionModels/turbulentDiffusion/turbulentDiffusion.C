/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
namespace mixingSubModels
{
namespace mixingDiffusionModels
{
    defineTypeNameAndDebug(turbulentDiffusion, 0);

    addToRunTimeSelectionTable
    (
        mixingDiffusionModel,
        turbulentDiffusion,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingDiffusionModels::turbulentDiffusion
::turbulentDiffusion
(
    const dictionary& dict
)
:
    mixingDiffusionModel(dict),
    gammaLam_(dict.lookup("gammaLam")),
    Sc_(readScalar(dict.lookup("Sc")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingDiffusionModels::turbulentDiffusion
::~turbulentDiffusion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::mixingSubModels::mixingDiffusionModels::turbulentDiffusion
::momentDiff
(
    const volScalarField& moment
) const
{
    volScalarField gamma(turbViscosity(moment)/Sc_ + gammaLam_);

    return fvm::laplacian(gamma, moment);
}

Foam::tmp<Foam::volScalarField>
Foam::mixingSubModels::mixingDiffusionModels::turbulentDiffusion::
turbViscosity(const volScalarField& moment) const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (moment.mesh().foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            moment.mesh().lookupObject<cmpTurbModel>
            (
                cmpTurbModel::propertiesName
            );

        return turb.mut()/turb.rho();
    }
    else if
    (
        moment.mesh().foundObject<icoTurbModel>(icoTurbModel::propertiesName)
    )
    {
        const incompressible::turbulenceModel& turb =
            moment.mesh().lookupObject<icoTurbModel>
            (
                icoTurbModel::propertiesName
            );

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
