/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
     \\/     M anipulation  |
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
    PrinceAndBlanch

Description
    Brownian aggregation kernel supplemented for particles in turbulent
    flows in conditions when the particle size is smaller than the local
    Kolmogorov scale.

    \f[
        \omega = C_a \epsilon_f^{1/3}(\xi_i^{2/3} + \xi_j^{2/3})
            \text{exp}
            \left(
              - C_{eff} \frac{\mu_f \rho_f \epsilon_f}{\sigma^2}
                \left(
                    \frac{\xi_i^{1/3} \xi_j^{1/3}}{\xi_i^{1/3} + \xi_j^{1/3}}
                \right)
            \right)
    \f]

    where

    \vartable
        \sigma      |  Surface tension
        \rho_f      |  Density of liquid phase [kg/m3]
        \epsilon_f  |  Turbulent dissipation rate of liquid phase [m2/s3]
        \mu_f       |  Molecular dynamic iscosity of liquid phase [Pa s]
        \xi_i       |  Abscissa of bubble i [m]
        \xi_j       |  Abscissa of bubble j [m]
    \endvartable

    References
    \verbatim
        "Description of interaction processes in agitated liquid-liquid dis[ersions"
        C Coulaloglou, L Tavlarides
        Chemical Engineering Science
        Volume 32, Pages 1289-1297, 1977
    \endverbatim

SourceFiles
    PrinceAndBlanch.C

\*---------------------------------------------------------------------------*/

#ifndef PrinceAndBlanch_H
#define PrinceAndBlanch_H

#include "coalesenceFrequencyKernel.H"
#include "twoPhaseSystem.H"
#include "phaseModel.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace coalesenceFrequencyKernels
{

/*---------------------------------------------------------------------------*\
                    Class PrinceAndBlanch Declaration
\*---------------------------------------------------------------------------*/

class PrinceAndBlanch
:
    public coalesenceFrequencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Switch for turbulence driven coalesence
        Switch turbulent_;

        //- Switch for buoyancy drive coalesence
        Switch buoyant_;

        //- Switch for laminar shear driven coalesence
        Switch LS_;


public:

        //- Runtime type information
        TypeName("PrinceAndBlanch");


    // Constructors

        //- Construct from components
        PrinceAndBlanch
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~PrinceAndBlanch();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> omega
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalesenceFrequencyKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
