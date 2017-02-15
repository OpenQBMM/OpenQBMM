/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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

#include "AyaziShamlou.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(AyaziShamlou, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        AyaziShamlou,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::AyaziShamlou
::AyaziShamlou
(
    const dictionary& dict
)
:
    breakupKernel(dict),
    A_(dict.lookup("A")),
    df_(dict.lookup("df")),
    H0_(dict.lookup("H0")),
    primarySize_(dict.lookup("primarySize"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::AyaziShamlou::~AyaziShamlou()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::breakupKernels::AyaziShamlou::Kb
(
    const volScalarField& abscissa
) const
{
    if (!abscissa.db().foundObject<fluidThermo>(basicThermo::dictName))
    {
        FatalErrorInFunction
            << "No valid thermophysical model found."
            << abort(FatalError);
    }

    const fluidThermo& flThermo =
        abscissa.db().lookupObject<fluidThermo>(basicThermo::dictName);

    typedef compressible::turbulenceModel cmpTurbModel;

    if
    (
        !abscissa.db().foundObject<cmpTurbModel>
        (
            cmpTurbModel::propertiesName
        )
    )
    {
        FatalErrorInFunction
            << "No valid compressible turbulence model found."
            << abort(FatalError);
    }

    const cmpTurbModel& flTurb =
        abscissa.db().lookupObject<cmpTurbModel>
        (
            turbulenceModel::propertiesName
        );

    // Interparticle force
    dimensionedScalar F = A_*primarySize_/(12.0*sqr(H0_));

    // Coefficient of volume fraction (Vanni, 2000)
    dimensionedScalar C = 0.41*df_ - 0.211;

    // Volume fraction of solid within aggregates
    volScalarField phiL(C*pow(abscissa/primarySize_, df_ - 3.0));

    // Coordination number
    volScalarField kc(15.0*pow(phiL, 1.2));

    // Aggregation strength
    volScalarField sigma(9.0*kc*phiL*F/(8.0*sqr(primarySize_)
            *Foam::constant::mathematical::pi));

    volScalarField epsilonByNu(flTurb.epsilon()/flThermo.nu());

    volScalarField tau(flThermo.mu()*sqrt(epsilonByNu));

    return sqrt(epsilonByNu/15.0)*exp(-sigma/tau);
}

// ************************************************************************* //
