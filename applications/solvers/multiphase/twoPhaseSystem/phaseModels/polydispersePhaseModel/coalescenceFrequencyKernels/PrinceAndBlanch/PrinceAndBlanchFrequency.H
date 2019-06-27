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
    Coalesence frequency kernel based on turbulence, buoyancy, and
    laminar shear.

    \f[
        \omega_{turb} =
            0.089 \pi (d_i + d_j)^2
            \sqrt{d_i^{\frac{2}{3}} + d_j^{\frac{2}{3}} }
            epsilon_f^{frac{1}{3}}

        \omega_{buoyancy} =
            \pi (d_i + d_j)^2
            (
               \sqrt{2.14*\sigma/(d_i \rho_f) + 0.5 \mathbf{g} d_i}
             - \sqrt{2.14*\sigma/(d_j \rho_f) + 0.5 \mathbf{g} d_j}
            )

        \omega_{shear} =
            \frac{2}{3} (d_i + d_j)^3 |\nabla \mathbf{U}_f|

    \f]

    where

    \vartable
        \sigma      |  Surface tension
        \mathbf{U}_f|  Liquid velocity [m/s]
        \rho_f      |  Density of liquid phase [kg/m3]
        \epsilon_f  |  Turbulent dissipation rate of liquid phase [m2/s3]
        g           |  Gravitational acceleration [m/s2]
        \d_i        |  Diameter of bubble i [m]
        \d_j        |  Diameter of bubble j [m]
    \endvartable

    References
    \verbatim
        "Bubble Coalescence and Break-up n Air-Sparged Bubble Columns"
        M. J. Prince, H. W. Blanch
        AIChE Journal
        Volume 36, Pages 1485 - 1499, 1990
    \endverbatim

SourceFiles
    PrinceAndBlanchFrequency.C

\*---------------------------------------------------------------------------*/

#ifndef PrinceAndBlanchFrequency_H
#define PrinceAndBlanchFrequency_H

#include "coalescenceFrequencyKernel.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace coalescenceFrequencyKernels
{

/*---------------------------------------------------------------------------*\
                    Class PrinceAndBlanch Declaration
\*---------------------------------------------------------------------------*/

class PrinceAndBlanch
:
    public coalescenceFrequencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Turbulent dissipation
        volScalarField epsilonf_;

        //- Switch for turbulence driven coalescence
        Switch turbulent_;

        //- Switch for buoyancy drive coalescence
        Switch buoyant_;

        //- Switch for laminar shear driven coalescence
        Switch LS_;


public:

        //- Runtime type information
        TypeName("PrinceAndBlanch");


    // Constructors

        //- Construct from components
        PrinceAndBlanch
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& continuousPhase
        );


    //- Destructor
    virtual ~PrinceAndBlanch();


    // Member Functions

        //- Update fields
        virtual void update
        (
            const fluidThermo& thermo,
            const turbulenceModel& turb
        );

        //- Aggregation kernel
        virtual scalar omega
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceFrequencyKernels
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
