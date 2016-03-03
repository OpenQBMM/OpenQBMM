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

#include "CoulaloglouTavlarides.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(CoulaloglouTavlarides, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        CoulaloglouTavlarides,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::CoulaloglouTavlarides
::CoulaloglouTavlarides
(
    const dictionary& dict
)
:
    breakupKernel(dict),
    C1b_
    (
        dict.lookupOrDefault
        (
            "C1b", 
            dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.00481)  
        )
    ),
    C2b_
    (
        dict.lookupOrDefault
        (
            "C2b", 
            dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.08)  
        )
    ),
    surfaceT_
    (
        dict.lookup("surfaceT")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::CoulaloglouTavlarides::~CoulaloglouTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::breakupKernels::CoulaloglouTavlarides::Kb
(
    const volScalarField& abscissa
) const
{   
    const compressible::turbulenceModel& fluidTurb =
        abscissa.mesh().lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

        const volScalarField& rhoc = fluidTurb.rho();


        return C1b_*pow(fluidTurb.epsilon(), 1.0/3.0)*pow(abscissa, -2.0/3.0)
            *exp(-C2b_*surfaceT_/rhoc/pow(fluidTurb.epsilon(), 2.0/3.0)/pow(abscissa, 5.0/3.0));
}

// ************************************************************************* //
