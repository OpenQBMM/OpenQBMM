/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
    Foam::coalescenceEfficiencyKernels::PrinceAndBlanch

Description
    Coalesence efficiency based on film drainage and rupture times.

    \f[
        P_c = \mathrm{exp}\left(\frac{t_{ij}}{\tau_{ij}}\right)
        t_{ij} = \sqrt{frac{\rho_f r_{ij}^3}{16 \sigma}} ln(\frac{h_0}{h_f})
        \tau_{ij} = \frac{r_{ij}^{2/3}}{\epsilon_f^{1/3}}
        r_{ij} = \frac{1}{2} \left(\frac{1}{r_i} + \frac{1}{r_j} \right)^{-1}
    \f]

    where

    \vartable
        \sigma      |  Surface tension
        \rho_f      |  Density of liquid phase [kg/m3]
        \epsilon_f  |  Turbulent dissipation rate of liquid phase [m2/s3]
        \h_f        |  Film thickness [m]
        \h_0        |  Initial film thickness [m]
        r_i         |  Radius of bubble i [m]
        r_j         |  Radius of bubble j [m]
    \endvartable

    References
    \verbatim
        "Bubble Coalescence and Break-up n Air-Sparged Bubble Columns"
         M Prince and Harvey Blanch
         M. J. Prince, H. W. Blanch
         AIChE Journal
         Volume 36, Pages 1485 - 1499, 1990
    \endverbatim

SourceFiles
    PrinceAndBlanchEfficiency.C

\*---------------------------------------------------------------------------*/

#ifndef PrinceAndBlanchEfficiency_H
#define PrinceAndBlanchEfficiency_H

#include "coalescenceEfficiencyKernel.H"
#include "twoPhaseSystem.H"
#include "phaseModel.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace coalescenceEfficiencyKernels
{

/*---------------------------------------------------------------------------*\
                    Class PrinceAndBlanch Declaration
\*---------------------------------------------------------------------------*/

class PrinceAndBlanch
:
    public coalescenceEfficiencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Initial film thickness
        dimensionedScalar ho_;

        //- Critical film thickness
        dimensionedScalar hf_;

        //- Turbulent dissipation
        volScalarField epsilonf_;


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
        virtual scalar Pc
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceEfficiencyKernels
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
