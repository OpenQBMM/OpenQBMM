/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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

#include "ChengVigilFoxFNP.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(ChengVigilFoxFNP, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        ChengVigilFoxFNP,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::ChengVigilFoxFNP
::ChengVigilFoxFNP
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    flThermo_(mesh_.lookupObject<fluidThermo>(basicThermo::dictName)),
    flTurb_
    (
        mesh_.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    ),
    T_(flThermo_.T()),
    rho_(flThermo_.rho()),
    mu_(flThermo_.mu()),
    epsilon_(flTurb_.epsilon()),
    mixingQuadrature_
    (
        mesh_.lookupObject<univariateQuadratureApproximation>
        (
            "mixingZZZ"
        )
    ),
    mixtureFraction_
    (
        mixingQuadrature_.moments()[1]
    ),
    molarMass_
    (
        dict.lookup("molarMass")
    ),
    molarVol1_
    (
        dict.lookup("molarVol1")
    ),
    molarVol2_
    (
        dict.lookup("molarVol2")
    ),
    soluteConc_
    (
        dict.lookup("soluteConc")
    ),
    aggregationEfficiency_
    (
        dict.lookupOrDefault("aggregationEfficiency", 1.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::ChengVigilFoxFNP
::~ChengVigilFoxFNP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::ChengVigilFoxFNP::Ka
(
    const scalar& abscissa1,
    const scalar& abscissa2,
    const label celli
) const
{
    // Calculate molar fraction
    scalar xa =
        (mixtureFraction_[celli]/molarVol1_.value())
       /(mixtureFraction_[celli]/molarVol1_.value()
      + (1.0 - mixtureFraction_[celli])/molarVol2_.value());

    // Local solute concentration, kmol/m3 -- InSoluteConc, mg/mL
    scalar cpcl =
        soluteConc_.value()*mixtureFraction_[celli]/molarMass_.value();

    // Equilibrium solute concentration, kmol/m3
    scalar ceq = 1200*exp(-14.533*(1 - xa))/molarMass_.value();

    // Flory parameters
    scalar FloryCoeff = 0.0064*exp(-3.15*xa);
    scalar FloryExp = 0.30 + 0.45*xa - 0.15*sqr(xa);

    // Radius of gyration of aggregating particle 1 (Flory's Law)
    scalar sqrt1
        = sqrt(FloryCoeff*pow(abscissa1*molarMass_.value(), 2.0*FloryExp));

    // Radius of gyration of aggregating particle 2 (Flory's Law)
    scalar sqrt2
        = sqrt(FloryCoeff*pow(abscissa2*molarMass_.value(), 2.0*FloryExp));

    // Brownian kernel
    scalar betaBrown =
         2.0*Foam::constant::physicoChemical::k.value()
        *T_[celli]*sqr(sqrt1 + sqrt2)/(3.0*mu_[celli]*sqrt1*sqrt2);

    //Turbulent kernel
    scalar betaTurb =
        1.2944*sqrt(epsilon_[celli]*rho_[celli]/mu_[celli])
       *pow3((1.0e-9)*(sqrt1 + sqrt2));

    return
        pos(cpcl - ceq)*Foam::constant::physicoChemical::NA.value()
       *(betaBrown + betaTurb);
}

// ************************************************************************* //
