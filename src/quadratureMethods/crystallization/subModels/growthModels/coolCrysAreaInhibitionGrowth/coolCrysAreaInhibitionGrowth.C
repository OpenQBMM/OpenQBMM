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

#include "coolCrysAreaInhibitionGrowth.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(coolCrysAreaInhibitionGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        coolCrysAreaInhibitionGrowth,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::growthModels::coolCrysAreaInhibitionGrowth
::coolCrysAreaInhibitionGrowth
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
    sigmaGrowthCrit_(dict.lookupOrDefault<scalar>("sigmaGrowthCrit", 0.0)),
    sigmaMax_(dict.lookupOrDefault<scalar>("sigmaMax", 10.0)),
    maxGrowth_
    (
        dict.lookupOrDefault
        (
            "maxGrowth",
            dimensionedScalar("maxGrowth", dimLength/dimTime, 1e-6)
        )
    ),
    m2Name_
    (
        dict.lookupOrDefault<word>("m2Name", "moment.2.populationBalance")
    ),
    alphaArea_
    (
        dict.lookupOrDefault
        (
            "alphaArea",
            dimensionedScalar("alphaArea", dimLength, 1e-4)
        )
    ),
    areaPower_(dict.lookupOrDefault<scalar>("areaPower", 1.0)),
    sigma_(mesh.lookupObject<volScalarField>("sigma")),
    activeCelli_(0),
    m2Ptr_(nullptr),
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
    if (!mesh.foundObject<volScalarField>(m2Name_))
    {
        FatalErrorInFunction
            << "Required second-moment field '" << m2Name_
            << "' was not found in the mesh object registry."
            << exit(FatalError);
    }

    m2Ptr_ = &mesh.lookupObject<volScalarField>(m2Name_);
}

Foam::populationBalanceSubModels::growthModels::coolCrysAreaInhibitionGrowth
::~coolCrysAreaInhibitionGrowth()
{}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysAreaInhibitionGrowth
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

    if (!std::isfinite(sigma_[celli]))
    {
        return 0.0;
    }

    const scalar sigmaEff =
        max(min(sigma_[celli], sigmaMax_) - sigmaGrowthCrit_, 0.0);

    if (sigmaEff <= SMALL)
    {
        return 0.0;
    }

    const scalar Gbase = Cg_.value()*Foam::pow(sigmaEff, powCoeff_);
    const scalar m2 = max((*m2Ptr_)[celli], 0.0);
    const scalar areaArg = max(alphaArea_.value()*m2, 0.0);

    scalar fArea =
        1.0/(1.0 + Foam::pow(areaArg, areaPower_));

    fArea = max(min(fArea, 1.0), 0.0);

    scalar G = Gbase*fArea;
    G = min(G, maxGrowth_.value());

    // This field is for visualization only. When Kg() is called for multiple
    // QMOM nodes in one cell, the stored value may correspond to the last
    // visited node rather than a node-weighted average.
    const_cast<volScalarField&>(GField_)[celli] = G;

    return G;
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysAreaInhibitionGrowth
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
Foam::populationBalanceSubModels::growthModels::coolCrysAreaInhibitionGrowth
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
