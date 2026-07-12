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

#include "coolCrysSizeEffectGrowth.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(coolCrysSizeEffectGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        coolCrysSizeEffectGrowth,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
::coolCrysSizeEffectGrowth
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
    Lc_
    (
        dict.lookupOrDefault
        (
            "Lc",
            dimensionedScalar("Lc", dimLength, 200e-6)
        )
    ),
    qSize_(dict.lookupOrDefault<scalar>("qSize", 2.0)),
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

Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
::~coolCrysSizeEffectGrowth()
{}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
::growthRateForLength
(
    const scalar& length,
    const label celli
) const
{
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

    const scalar Li = max(length, SMALL);
    const scalar Gbase = Cg_.value()*Foam::pow(sigmaEff, powCoeff_);

    scalar fSize =
        1.0/(1.0 + Foam::pow(Li/Lc_.value(), qSize_));

    fSize = max(min(fSize, 1.0), 0.0);

    scalar G = Gbase*fSize;
    G = min(G, maxGrowth_.value());

    return G;
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
::Kg
(
    const scalar& abscissa,
    const bool lengthBased,
    const label environment
) const
{
    const label celli = activeCelli_;
    (void)environment;

    const scalar Li =
        lengthBased
      ? abscissa
      :
        // When lengthBased is false, assume the abscissa is a volume-like
        // coordinate and map it to a characteristic length through cbrt.
        Foam::cbrt(max(abscissa, SMALL));

    const scalar G = growthRateForLength(Li, celli);
    const_cast<volScalarField&>(GField_)[celli] = G;

    return G;
}


void Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
::updateMeanGrowthRate
(
    const label celli,
    const scalarQuadratureApproximation& quadrature
) const
{
    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    const label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        const_cast<volScalarField&>(GField_)[celli] = 0.0;
        return;
    }

    scalar weightedGrowth = 0.0;
    scalar totalNumberDensity = 0.0;

    forAll(nodes, nodeI)
    {
        const volScalarNode& node = nodes[nodeI];
        const scalar bAbscissa =
            max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));
        const scalar d = node.d(celli, bAbscissa);
        const scalar numberDensity =
            node.n(celli, node.primaryWeight()[celli], bAbscissa);
        const scalar G = Kg(d, node.lengthBased());

        if (std::isfinite(numberDensity) && std::isfinite(G) && numberDensity > SMALL)
        {
            weightedGrowth += numberDensity*G;
            totalNumberDensity += numberDensity;
        }
    }

    const_cast<volScalarField&>(GField_)[celli] =
        totalNumberDensity > SMALL
      ? weightedGrowth/totalNumberDensity
      : scalar(0.0);
}


void Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
::updateMeanGrowthRate
(
    const label celli,
    const velocityQuadratureApproximation& quadrature
) const
{
    const mappedPtrList<volVelocityNode>& nodes = quadrature.nodes();
    const label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        const_cast<volScalarField&>(GField_)[celli] = 0.0;
        return;
    }

    scalar weightedGrowth = 0.0;
    scalar totalNumberDensity = 0.0;

    forAll(nodes, nodeI)
    {
        const volVelocityNode& node = nodes[nodeI];
        const scalar bAbscissa =
            max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));
        const scalar d = node.d(celli, bAbscissa);
        const scalar numberDensity =
            node.n(celli, node.primaryWeight()[celli], bAbscissa);
        const scalar G = Kg(d, node.lengthBased());

        if (std::isfinite(numberDensity) && std::isfinite(G) && numberDensity > SMALL)
        {
            weightedGrowth += numberDensity*G;
            totalNumberDensity += numberDensity;
        }
    }

    const_cast<volScalarField&>(GField_)[celli] =
        totalNumberDensity > SMALL
      ? weightedGrowth/totalNumberDensity
      : scalar(0.0);
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
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
Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth
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
