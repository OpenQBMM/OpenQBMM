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
    CoulaloglouAndTavlarides

Description
    Coalesence frequency kernel based on kinetic theory of gasses

    \f[
        \omega = \epsilon_f^{1/3}(d_i^{2/3} + d_j^{2/3})
    \f]

    where

    \vartable
        \epsilon_f  |  Turbulent dissipation rate of liquid phase [m2/s3]
        d_i         |  Diameter of bubble i [m]
        d_j         |  Diameter of bubble j [m]
    \endvartable

    References
    \verbatim
        "Description of interaction processes in agitated liquid-liquid dispersions"
        C Coulaloglou, L Tavlarides
        Chemical Engineering Science
        Volume 32, Pages 1289-1297, 1977
    \endverbatim

SourceFiles
    CoulaloglouAndTavlaridesFrequency.C

\*---------------------------------------------------------------------------*/

#ifndef CoulaloglouAndTavlaridesFrequency_H
#define CoulaloglouAndTavlaridesFrequency_H

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
                    Class CoulaloglouAndTavlarides Declaration
\*---------------------------------------------------------------------------*/

class CoulaloglouAndTavlarides
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
        TypeName("CoulaloglouAndTavlarides");


    // Constructors

        //- Construct from components
        CoulaloglouAndTavlarides
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& continuousPhase
        );


    //- Destructor
    virtual ~CoulaloglouAndTavlarides();


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
