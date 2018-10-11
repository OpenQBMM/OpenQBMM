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

#include "coulaloglouTavlarides.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace coalescenceKernels
{
    defineTypeNameAndDebug(coulaloglouTavlarides, 0);

    addToRunTimeSelectionTable
    (
        coalescenceKernel,
        coulaloglouTavlarides,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::coalescenceKernels::coulaloglouTavlarides
::coulaloglouTavlarides
(
    const dictionary& dict
)
:
    etaC_(dict.lookupOrDefault("etaC", 6.0E9)),
    coalescenceKernel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::coalescenceKernels::coulaloglouTavlarides
::~coulaloglouTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::coalescenceKernels::coulaloglouTavlarides::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{
    if (!abscissa1.mesh().foundObject<fluidThermo>(basicThermo::dictName))
    {
        FatalErrorInFunction
            << "No valid thermophysical model found."
            << abort(FatalError);
    }

    const fluidThermo& flThermo =
        abscissa1.mesh().lookupObject<fluidThermo>(basicThermo::dictName);

    typedef compressible::turbulenceModel cmpTurbModel;

    if
    (
        !abscissa1.mesh().foundObject<cmpTurbModel>
        (
            cmpTurbModel::propertiesName
        )
    )
    {
        FatalErrorInFunction
            << "No valid compressible turbulence model found."
            << abort(FatalError);
    }

    const compressible::turbulenceModel& flTurb =
        abscissa1.mesh().lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

                    // NEED TO LOOK UP SIGMA STILL\\
        
    dimensionedScalar smallAbs("smallAbs", sqr(abscissa1.dimensions()), SMALL);

    volScalarField eta = exp(-etaC_*sqrt(2.0*flThermo.mu()/flTurb.nu()
       *pow(flTurb.epsilon(), 2.0/3.0)/sigma*(abscissa1*abscissa2)
       /(max(smallAbs, abscissa1 + abscissa2))));
    
    return
        Cc_*pow(flTurb.epsilon(), 1.0/3.0)*sqr(abscissa1 + abscissa2)
       *sqrt(pow(abscissa1, 2.0/3.0) + pow(abscissa2, 2.0/3.0))*eta;
}

// ************************************************************************* //
