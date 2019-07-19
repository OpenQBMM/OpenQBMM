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
    Foam::populationBalanceSubModels::coalescenceEfficiencyKernels::
    CoulaloglouAndTavlarides

Description
    Coalesence efficiency based on film drainage and rupture times.

    \f[
        P_c =
            \text{exp}
            \left(
              - C_{eff} \frac{\mu_f \rho_f \epsilon_f}{\sigma^2}
                \left(
                    \frac{d_i^{1/3} d_j^{1/3}}{d_i^{1/3} + d_j^{1/3}}
                \right)
            \right)
    \f]

    where

    \vartable
        \sigma      |  Surface tension
        \rho_f      |  Density of liquid phase [kg/m3]
        \epsilon_f  |  Turbulent dissipation rate of liquid phase [m2/s3]
        \mu_f       |  Molecular dynamic iscosity of liquid phase [Pa s]
        d_i         |  Diameter of bubble i [m]
        d_j         |  Diameter of bubble j [m]
    \endvartable

    References
    \verbatim
        "Description of interaction processes in agitated liquid-liquid
        dispersions"
        C Coulaloglou, L Tavlarides
        Chemical Engineering Science
        Volume 32, Pages 1289-1297, 1977
    \endverbatim

SourceFiles
    CoulaloglouAndTavlaridesEfficiency.C

\*---------------------------------------------------------------------------*/

#ifndef CoulaloglouAndTavlaridesEfficiency_H
#define CoulaloglouAndTavlaridesEfficiency_H

#include "coalescenceEfficiencyKernel.H"
#include "twoPhaseSystem.H"

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
                    Class CoulaloglouAndTavlarides Declaration
\*---------------------------------------------------------------------------*/

class CoulaloglouAndTavlarides
:
    public coalescenceEfficiencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Efficiency coefficient
        dimensionedScalar Ceff_;

        //- Turbulent dissipation
        volScalarField epsilonf_;

        //- Kinematic viscosity of fluid phase
        volScalarField muf_;


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
