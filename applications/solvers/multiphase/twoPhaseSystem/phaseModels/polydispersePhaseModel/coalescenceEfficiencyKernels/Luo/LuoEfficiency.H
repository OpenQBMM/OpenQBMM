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
        P_c = \mathrm{exp}
            \left(
               -\frac{\frac{3}{4}(1+\xi^2)^{1/2}(1+\xi^3)^{1/2}}
                {\left(\frac{\rho_b}{\rho_f}+C_{vm}\right)^{1/2}(1+\xi)^3}
                \mathrm{We}^{1/2}
            \right)
        u_{rel} =  2 \epsilon_f^{1/3}(d_i^{2/3} + d_j^{2/3})^{1/2}
        \mathrm{We} = \frac{\mathrm{min}(d_i, d_j)\rho_f\u_{rel}^2}{\sigma}
    \f]

    where

    \vartable
        \sigma      |  Surface tension
        \rho_f      |  Density of liquid phase [kg/m3]
        \rho_b      |  Density of bubble phase [kg/m3]
        \epsilon_f  |  Turbulent dissipation rate of liquid phase [m2/s3]
        d_i         |  Diameter of bubble i [m]
        d_j         |  Diameter of bubble j [m]
    \endvartable

    References
    \verbatim
        "Coalescence, breakup and liquid circulation in bubble column reactors"
         Luo, H.,
         PhD Theses,
         The Norwegian Institute of Technology, 1993
    \endverbatim

SourceFiles
    LuoEfficiency.C

\*---------------------------------------------------------------------------*/

#ifndef LuoEfficiency_H
#define LuoEfficiency_H

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
                            Class Luo Declaration
\*---------------------------------------------------------------------------*/

class Luo
:
    public coalescenceEfficiencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Turbulent dissipation
        volScalarField epsilonf_;

        //- Virutal mass coefficient
        volScalarField Cvm_;


public:

        //- Runtime type information
        TypeName("Luo");


    // Constructors

        //- Construct from components
        Luo
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& continuousPhase
        );


    //- Destructor
    virtual ~Luo();


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
