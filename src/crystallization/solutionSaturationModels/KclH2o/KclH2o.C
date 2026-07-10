/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "KclH2o.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solutionSaturationModels
{
    defineTypeNameAndDebug(KclH2o, 0);
    addToRunTimeSelectionTable(solutionSaturationModel, KclH2o, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::KclH2o::KclH2o
(
    const dictionary& dict,
    const objectRegistry& db
)
:
    solutionSaturationModel(db)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::KclH2o::~KclH2o()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::solutionSaturationModels::KclH2o::Csat
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tCsat
    (
        volScalarField::New
        (
            "Csat",
            T.mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );

    volScalarField& Csat = tCsat.ref();

    // KCl solubility correlation: Csat = a + b*T + c*T^2
    // where T is in Celsius, result is in g KCl/100g H2O
    // Converting to kg/m³ (mass concentration in solution)
    // For KCl solution: C[kg/m³] = S * ρ_solution / (100 + S)
    // where S is solubility in g/100g water
    // Approximation for dilute solutions: C[kg/m³] ≈ 10 * S[g/100g water]
    
    const scalar rhoSolution = 1000.0;  // Solution density [kg/m³], can be refined
    
    forAll(Csat, celli)
    {
        scalar Tc = T[celli] - 273.15;  // Convert K to °C
        // Solubility in g/100g water: 27.983 + 0.3193*T - 0.0004*T^2
        scalar S = 27.983 + 0.3193*Tc - 0.0004*pow(Tc, 2);
        // Convert to kg/m³: C = (S/1000) * rho_solution * 100/(100+S)
        // Simplified: C ≈ S * rho_solution / 100 for S << 100
        Csat[celli] = (S * rhoSolution) / 100.0;
    }

    volScalarField::Boundary& CsatBf = Csat.boundaryFieldRef();

    forAll(Csat.boundaryField(), patchi)
    {
        scalarField& Csatp = CsatBf[patchi];
        const scalarField& Tp = T.boundaryField()[patchi];

        forAll(Csatp, facei)
        {
            scalar Tc = Tp[facei] - 273.15;  // Convert K to °C
            scalar S = 27.983 + 0.3193*Tc - 0.0004*pow(Tc, 2);
            Csatp[facei] = (S * rhoSolution) / 100.0;
        }
    }

    return tCsat;
}


// ************************************************************************* //
