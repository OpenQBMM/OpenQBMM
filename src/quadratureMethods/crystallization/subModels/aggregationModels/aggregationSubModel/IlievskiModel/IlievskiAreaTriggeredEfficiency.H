/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Junjie Li
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

Class
    Foam::populationBalanceSubModels::aggregationKernels::
    crystalAggregationEfficiencies::IlievskiAreaTriggered

Description
    Ilievski aggregation efficiency with an additive area-triggered boost.

    The base efficiency is the existing Ilievski model:
    \f[
        \Psi_0 = \frac{0.00042}{1 + 6.25 \times 10^{-19} \mu^{1.5}
            \left( \frac{(\gamma^2 / U_{tip})^{2.25}}{(G/L_{10})^3} \right)}
    \f]

    The area-triggered activation adds
    \f[
        activation = 1 - exp(-(alphaAreaAgg * max(m_2 - m2Offset, 0))^{areaPowerAgg})
    \f]
    and
    \f[
        \Psi = min(1, \Psi_0 + psiBoostMax * activation)
    \f]

SourceFiles
    IlievskiAreaTriggeredEfficiency.C

\*---------------------------------------------------------------------------*/

#ifndef IlievskiAreaTriggeredEfficiency_H
#define IlievskiAreaTriggeredEfficiency_H

#include "crystalAggregationEfficiency.H"
#include "turbulentTransportModel.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalAggregationEfficiencies
{

class IlievskiAreaTriggered
:
    public crystalAggregationEfficiency
{
    // Private data

        scalar Utip_;
        scalar L10_;
        bool useConfiguredL10_;
        word m2Name_;
        scalar m2Offset_;
        scalar alphaAreaAgg_;
        scalar areaPowerAgg_;
        scalar psiBoostMax_;

        const volScalarField& GField_;
        const volScalarField& nu_;
        const volScalarField& epsilon_;
        const volScalarField& m2Field_;


public:

    TypeName("IlievskiAreaTriggered");

    IlievskiAreaTriggered
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~IlievskiAreaTriggered();

    virtual scalar Pc
    (
        const scalar& d1,
        const scalar& d2,
        const vector& Ur,
        const label celli
    ) const;
};

} // End namespace crystalAggregationEfficiencies
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

#endif

// ************************************************************************* //
