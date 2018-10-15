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
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Luo();


    // Member Functions

        //- Update fields
        virtual void update();

        //- Aggregation kernel
        virtual tmp<volScalarField> Pc
        (
            const label nodei,
            const label nodej
        ) const;

        //- Aggregation kernel
        virtual scalar Pc
        (
            const label celli,
            const label nodei,
            const label nodej
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceEfficiencyKernels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
