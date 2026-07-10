/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 Alberto Passalacqua
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

#include "coolCrysGrowth.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(coolCrysGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        coolCrysGrowth,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::growthModels::coolCrysGrowth
::coolCrysGrowth
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    growthModel(dict, mesh),
    Cg_
    (
        dict.lookupOrDefault
        (
            "CgGrowth",
            dimensionedScalar("CgGrowth", dimLength/dimTime, 1e-6)
        )
    ),
    powCoeff_(readScalar(dict.lookup("powCoeff"))),
    sigma_(mesh.lookupObject<volScalarField>("sigma")),
    GField_
    (
        IOobject
        (
            "growthRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength/dimTime, 0.0)
    )
{
}

Foam::populationBalanceSubModels::growthModels::coolCrysGrowth
::~coolCrysGrowth()
{}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysGrowth
::Kg
(
    const scalar& abscissa,
    const bool lengthBased,
    const label celli,
    const label environment
) const
{

    scalar G = 0.0;
    
    // ENHANCED PROTECTION: Check for numerical stability
    // Reject unrealistic sigma values that indicate numerical problems
    const scalar sigmaMax = 10.0;
    const scalar sigmaMin = 1e-10;
    
    if (sigma_[celli] > sigmaMin && sigma_[celli] < sigmaMax)
    {
        G = Cg_.value() * pow(sigma_[celli], powCoeff_);
        const scalar maxG = 1.0e-2;
        G = min(G, maxG);
    }
    else if (sigma_[celli] >= sigmaMax)
    {
        G = Cg_.value() * pow(sigmaMax, powCoeff_);
        G = min(G, 1.0e-3);
    }
    // else: sigma too small or negative -> G = 0 (no growth)
    
    // Update GField for output visualization
    const_cast<volScalarField&>(GField_)[celli] = G;


    return G;          // dL/dt

}

// ************************************************************************* //
