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
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(CoulaloglouTavlarides, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        CoulaloglouTavlarides,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::CoulaloglouTavlarides::CoulaloglouTavlarides
(
    const dictionary& dict
)
:
    aggregationKernel(dict),
    Ca2_
    (
        dict.lookupOrDefault
        (
            "Ca2", 
            dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.88)  
        )
    ),
    surfaceT_
    (
        dict.lookup("surfaceT")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::CoulaloglouTavlarides
::~CoulaloglouTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::aggregationKernels::CoulaloglouTavlarides::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{   
    if (!abscissa1.mesh().foundObject<fluidThermo>(basicThermo::dictName))
    {
        FatalErrorIn
        (
            "Foam::populationBalanceSubModels::aggregationKernels::"
            "CoulaloglouTavlarides::Ka\n"
            "(\n"
            "   const volScalarField& abscissa1,\n"
            "   const volScalarField& abscissa2,\n"
            ")"
        )   << "No valid thermophysical model found."
            << abort(FatalError);
    }
    
    typedef compressible::turbulenceModel cmpTurbModel;
    
    if 
    (
        !abscissa1.mesh().foundObject<cmpTurbModel>
        (
            cmpTurbModel::propertiesName
        )
    )
    {
        FatalErrorIn
        (
            "Foam::populationBalanceSubModels::aggregationKernels::"
            "CoulaloglouTavlarides::Ka\n"
            "(\n"
            "   const volScalarField& abscissa1,\n"
            "   const volScalarField& abscissa2,\n"
            ")"
        )   << "No valid compressible turbulence model found."
            << abort(FatalError);
    }
        
    const compressible::turbulenceModel& flTurb =
        abscissa1.mesh().lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

        dimensionedScalar D("D", dimensionSet(0, -2, 0, 0, 0), 1.0);

        return Ca2_*Foam::pow(flTurb.epsilon(),1.0/3.0) * Foam::pow(abscissa1 + abscissa2, 2.0)
        * Foam::pow(Foam::pow(abscissa1, 2.0/3.0) 
            + Foam::pow(abscissa2, 2.0/3.0), 1.0/2.0)
        * Foam::exp(-Ca_*D*flTurb.mu()*flTurb.rho()*flTurb.epsilon()/Foam::pow(surfaceT_, 2.0)
            *Foam::pow(abscissa2*abscissa1/(abscissa1+abscissa2), 4.0));
}

// ************************************************************************* //
