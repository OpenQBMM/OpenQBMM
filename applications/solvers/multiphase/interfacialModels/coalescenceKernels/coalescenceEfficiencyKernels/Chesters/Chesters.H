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
    Foam::populationBalanceSubModels::coalescenceEfficiencyKernels::
    Chesters

Description
    Coaelescence efficieny kernel based on the interaction time of bubble.
    Only valid for bubbles of equal sizes.

    \f[
        P_c =
            \text{exp}
            (
            - \frac{\theta \sqrt{We} \sqrt{0.75 (1.0 + \xi^2) (1.0 + \xi^3)}}
              {\sqrt{\rho_b/rho_f}*(1.0 + \xi)^3}
            );
        \xi = \frac{d_i}{d_j}
        \theta = x Re^{y} \text{We}^{z}
    \f]

    where

    \vartable
        \rho_f      |  Density of liquid phase [kg/m3]
        \rho_b      |  Density of bubble phase [kg/m3]
        We          |  Webber number [ ]
        Re          |  Reynolds number [ ]
        x           |  Efficiency cofficient [ ]
        y           |  Reynolds number coefficient [ ]
        z           |  Weber number coefficient [ ]
        d_i         |  Diameter of bubble i [m]
        d_j         |  Diamerer of bubble j [m]
    \endvartable

    References
    \verbatim
        "The Modelling of coalescence processes in fluid-liquid dispersions:
        a review of current understanding"
        A.K. Chesters
        Chemical Engineering Research and Design
        Volume 69, Pages 259 - 270, 1991
    \endverbatim

SourceFiles
    Chesters.C

\*---------------------------------------------------------------------------*/

#ifndef Chesters_H
#define Chesters_H

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
                    Class Chesters Declaration
\*---------------------------------------------------------------------------*/

class Chesters
:
    public coalescenceEfficiencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Efficiency coefficient
        dimensionedScalar Ceff_;

        //- Reynolds number exponent
        dimensionedScalar ReExp_;

        //- Webber number exponent
        dimensionedScalar WeExp_;


public:

        //- Runtime type information
        TypeName("Chesters");


    // Constructors

        //- Construct from components
        Chesters
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Chesters();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> Pc
        (
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
