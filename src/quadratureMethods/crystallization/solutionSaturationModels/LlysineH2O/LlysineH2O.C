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

#include "LlysineH2O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solutionSaturationModels
{
    defineTypeNameAndDebug(LlysineH2O, 0);
    addToRunTimeSelectionTable(solutionSaturationModel, LlysineH2O, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::LlysineH2O::LlysineH2O
(
    const dictionary& dict,
    const objectRegistry& db
)
:
    solutionSaturationModel(db)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::LlysineH2O::~LlysineH2O()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::solutionSaturationModels::LlysineH2O::Csat
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

    // Coefficients for L-Lysine Free Base Solubility
    // Fit Model: S = A + B*T + C*T^2
    // Unit: kg/m^3
    // Temperature: Kelvin
    const scalar A = 3041.7992;
    const scalar B = -30.1725;
    const scalar C = 0.0750;

    // 1. Calculate Internal Field
    forAll(Csat, celli)
    {
        // T is already in Kelvin in OpenFOAM
        scalar TK = T[celli]; 
        
        // Direct calculation in kg/m^3
        // S = 3041.80 - 30.173*T + 0.075*T^2
        Csat[celli] = A + B*TK + C*pow(TK, 2);
    }

    // 2. Calculate Boundary Field
    volScalarField::Boundary& CsatBf = Csat.boundaryFieldRef();

    forAll(Csat.boundaryField(), patchi)
    {
        scalarField& Csatp = CsatBf[patchi];
        const scalarField& Tp = T.boundaryField()[patchi];

        forAll(Csatp, facei)
        {
            scalar TK = Tp[facei];
            
            // Direct calculation in kg/m^3
            Csatp[facei] = A + B*TK + C*pow(TK, 2);
        }
    }

    return tCsat;
}


// ************************************************************************* //