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
    Foam::populationBalanceSubModels::aggregationKernels::linnikovAggregation

Description
    Linnikov aggregation kernel for cooling crystallization based on 
    supersaturation-dependent aggregation.

    The kernel assumes size-independent aggregation and is calculated as:
    \f[
        K_a = K_1 \left( 1 - e^{-D_3 \cdot e^{-\frac{D_1}{\sigma + D_2}}} \right)
    \f]

    Or in simplified form (first-order expansion):
    \f[
        K_a = K_1 D_3 e^{-\frac{D_1}{\sigma + D_2}}
    \f]

    where \f$\sigma = \ln(C/C_0)\f$ is the relative supersaturation.

    Model parameters:
    \f[
        D_1 = \frac{4 h E_\sigma^2 N_A M}{\rho_c R^2 T^2}
    \f]
    \f[
        D_2 = \frac{2 E_\sigma M}{h \rho_c R T}
    \f]
    \f[
        D_3 = I_0 t_k
    \f]

    \vartable
        h           | Height of nucleus-bridge [m]
        M           | Molar weight of crystal [kg/mol]
        \rho_c      | Crystal density [kg/m^3]
        R           | Universal gas constant [J/(mol·K)]
        T           | Temperature [K]
        N_A         | Avogadro's number [1/mol]
        E_\sigma    | Specific surface energy [J/m^2]
        I_0         | Rate of nucleus-bridge formation [1/s]
        t_k         | Contact time [s]
    \endvartable

    Reference values from Linnikov (2008):
    - K2SO4: D1 = 0.1459, D2 = 0.044, K1*D3 = 9.5446e-8
    - KNO3:  D1 = 0.1459, D2 = 0.044, K1*D3 = 4.4544e-7

    Reference:
    \verbatim
        Linnikov, O.D. (2008).
        "Mechanism of aggregation and intergrowth of crystals during 
        bulk crystallization from solutions"
        Crystal Research and Technology
    \endverbatim

SourceFiles
    linnikovAggregation.C

\*---------------------------------------------------------------------------*/

#ifndef linnikovAggregation_H
#define linnikovAggregation_H

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class linnikovAggregation Declaration
\*---------------------------------------------------------------------------*/

class linnikovAggregation
:
    public aggregationKernel
{
    // Private data

        //- Parameter D1: related to surface energy and temperature
        //  D1 = 4*h*E_sigma^2*N_A*M / (rho_c*R^2*T^2)
        scalar D1_;

        //- Parameter D2: related to surface energy
        //  D2 = 2*E_sigma*M / (h*rho_c*R*T)
        scalar D2_;

        //- Combined parameter K1*D3
        //  K1*D3 = K1 * I0 * tk
        scalar K1D3_;

        //- Use simplified (first-order) model
        bool simplified_;

        //- Reference to the solver supersaturation field
        const volScalarField& sigma_;


public:

    //- Runtime type information
    TypeName("linnikovAggregation");


    // Constructors

        //- Construct from components
        linnikovAggregation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~linnikovAggregation();


    // Member Functions

        //- Aggregation kernel (size-independent)
        virtual scalar Ka
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli,
            const label environment = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
