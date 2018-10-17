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
    Foam::populationBalanceSubModels::aggregationKernels::Luo

Description
    Coalesence frequency kernel based on turbulence, buoyancy, and
    laminar shear.

    \f[
        \omega = \frac{\pi}{4}(d_i + d_j)^2 u_{rel}

        u_{rel} =  2 \epsilon_f^{1/3}(d_i^{2/3} + d_j^{2/3})^{1/2}
    \f]

    where

    \vartable
        \epsilon_f  |  Turbulent dissipation rate of liquid phase [m2/s3]
        \d_i        |  Diameter of bubble i [m]
        \d_j        |  Diameter of bubble j [m]
    \endvartable

    References
    \verbatim
        "Coalescence, breakup and liquid circulation in bubble column reactors"
         Luo, H.,
         PhD Theses,
         The Norwegian Institute of Technology, 1993
    \endverbatim

SourceFiles
    LuoFrequecy.C

\*---------------------------------------------------------------------------*/

#ifndef PrinceAndBlanch_H
#define PrinceAndBlanch_H

#include "coalescenceFrequencyKernel.H"
#include "twoPhaseSystem.H"
#include "phaseModel.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceFrequencyKernels
{

/*---------------------------------------------------------------------------*\
                            Class Luo Declaration
\*---------------------------------------------------------------------------*/

class Luo
:
    public coalescenceFrequencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Turbulent dissipation
        volScalarField epsilonf_;


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
        virtual tmp<volScalarField> omega
        (
            const label nodei,
            const label nodej
        ) const;

        //- Aggregation kernel
        virtual scalar omega
        (
            const label celli,
            const label nodei,
            const label nodej
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceFrequencyKernels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
