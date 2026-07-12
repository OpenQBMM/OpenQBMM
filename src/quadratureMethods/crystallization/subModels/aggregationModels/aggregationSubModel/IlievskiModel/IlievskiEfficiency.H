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
    crystalAggregationEfficiencies::Ilievski

Description
    Ilievski aggregation efficiency model for cooling crystallization.

    The efficiency is calculated as:
    \f[
        \Psi_{10} = \frac{0.00042}{1 + 6.25 \times 10^{-19} \mu^{1.5}
            \left( \frac{(\gamma^2 / U_{tip})^{2.25}}{(G/L_{10})^3} \right)}
    \f]

    \vartable
        \Psi_{10}   | Number-based aggregation efficiency
        \mu         | Dynamic viscosity [Pa.s]
        \gamma      | Shear rate [1/s]
        U_{tip}     | Impeller tip speed [m/s]
        G           | Crystal growth rate [m/s]
        L_{10}      | Number-based mean crystal size [m]
    \endvartable

    Reference:
    \verbatim
        Ilievski, D. (2006).
        "Agglomeration Efficiency Model"
        DOI: ...
    \endverbatim

SourceFiles
    IlievskiEfficiency.C

\*---------------------------------------------------------------------------*/

#ifndef IlievskiEfficiency_H
#define IlievskiEfficiency_H

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
                         Class Ilievski Declaration
\*---------------------------------------------------------------------------*/

class Ilievski
:
    public crystalAggregationEfficiency
{
    // Private data

        scalar Utip_;
        scalar psi0_;
        scalar CIlievski_;
        scalar rho_;
        scalar L10_;
        bool useConfiguredL10_;

        const volScalarField& GField_;
        const volScalarField& nu_;
        const volScalarField& epsilon_;


public:

    //- Runtime type information
    TypeName("Ilievski");


    // Constructors

        //- Construct from components
        Ilievski
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Ilievski();


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
