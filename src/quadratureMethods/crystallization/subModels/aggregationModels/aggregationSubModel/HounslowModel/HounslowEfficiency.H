/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Junjie Li
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
    crystalAggregationEfficiencies::Hounslow

Description
    Hounslow aggregation efficiency model for cooling crystallization.

    The efficiency is calculated as:
    \f[
        \psi(M) = \frac{1}{2}\left[1 + \mathrm{erf}\left(
            \frac{\ln(M/\bar{M}_g)}{\sqrt{2}\ln\sigma_g}
        \right)\right]
    \f]

    where
    \f[
        M = \frac{\sigma_Y G L^*}{\dot{\gamma}^2 \mu L^2}
    \f]

    \vartable
        M           | Dimensionless parameter
        \bar{M}_g   | Geometric mean of M
        \sigma_g    | Geometric standard deviation
        \sigma_Y    | Yield stress
        G           | Crystal growth rate [m/s]
        L^*         | Characteristic length scale
        \dot{\gamma}| Shear rate [1/s]
        \mu         | Dynamic viscosity [Pa.s]
        L           | Crystal size [m]
    \endvartable

    Reference:
    \verbatim
        Hounslow, M.J. et al. (2013).
        "Aggregation During Growing Crystals"
        DOI: ...
    \endverbatim

SourceFiles
    HounslowEfficiency.C

\*---------------------------------------------------------------------------*/

#ifndef HounslowEfficiency_H
#define HounslowEfficiency_H

#include "crystalAggregationEfficiency.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalAggregationEfficiencies
{

/*---------------------------------------------------------------------------*\
                         Class Hounslow Declaration
\*---------------------------------------------------------------------------*/

class Hounslow
:
    public crystalAggregationEfficiency
{
    // Private data

        //- Geometric mean of M
        scalar Mg_;

        //- Geometric standard deviation
        scalar sigmag_;

        //- Composite material parameter L*sigma_Y [N/m]
        //  Typical values: Rounded COM: 1.35, Well-defined COM: 3.68
        //  Calcite with KNO3: 5.05, Calcite without KNO3: 1.54
        scalar Lstar_sigmaY_;

        //- Fluid density used to convert kinematic to dynamic viscosity
        scalar rhoFluid_;

        //- q-width over which Hounslow contact efficiency is smoothly tapered
        scalar qBlend_;

        //- Lower bound used only for contact-length evaluation inside the taper
        scalar qMin_;

        //- Reference to growth rate field (from growthModel)
        const volScalarField& GField_;

        //- Reference to kinematic viscosity field
        const volScalarField& nu_;

        //- Reference to turbulent dissipation rate
        const volScalarField& epsilon_;


public:

    //- Runtime type information
    TypeName("Hounslow");


    // Constructors

        //- Construct from components
        Hounslow
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Hounslow();


    // Member Functions

        //- Calculate aggregation efficiency
        virtual scalar Pc
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace crystalAggregationEfficiencies
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
