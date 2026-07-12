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

#include "coolCrysGrowthDissolution.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(coolCrysGrowthDissolution, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        coolCrysGrowthDissolution,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::growthModels::coolCrysGrowthDissolution
::coolCrysGrowthDissolution
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
    CgDiss_
    (
        dict.lookupOrDefault
        (
            "CgDiss",
            dimensionedScalar("CgDiss", dimLength/dimTime, 1e-6)
        )
    ),
    powCoeffDiss_(readScalar(dict.lookup("powCoeffDiss"))),
    sigma_(mesh.lookupObject<volScalarField>("sigma")),
    activeCelli_(0),
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

Foam::populationBalanceSubModels::growthModels::coolCrysGrowthDissolution
::~coolCrysGrowthDissolution()
{}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysGrowthDissolution
::Kg
(
    const scalar& abscissa,
    const bool lengthBased,
    const label environment
) const
{
    const label celli = activeCelli_;
    (void)abscissa;
    (void)lengthBased;
    (void)environment;

    scalar G = 0.0;

    const scalar limitRate = 1.0e-2; 
    const scalar sigmaMax = 10.0;
    const scalar sigmaMin = 1e-10;
    
    scalar sigmaVal = sigma_[celli];

    if (sigmaVal > sigmaMin)
    {
        scalar effSigma = min(sigmaVal, sigmaMax);
        
        G = Cg_.value() * pow(effSigma, powCoeff_);
        
        G = min(G, limitRate);
    }
    else if (sigmaVal < -sigmaMin)
    {
        scalar absSigma = mag(sigmaVal);
    
        scalar effSigma = min(absSigma, sigmaMax);
        
        G = -1.0 * CgDiss_.value() * pow(effSigma, powCoeffDiss_);
        
        G = max(G, -limitRate);
    }
    else
    {
        G = 0.0;
    }
    const_cast<volScalarField&>(GField_)[celli] = G;
    
    return G;
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysGrowthDissolution
::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    activeCelli_ = celli;
    return growthModel::phaseSpaceConvection(momentOrder, celli, quadrature);
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysGrowthDissolution
::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    activeCelli_ = celli;
    return growthModel::phaseSpaceConvection(momentOrder, celli, quadrature);
}
// ************************************************************************* //
