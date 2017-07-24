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
    Foam::populationBalanceSubModels::coalesenceEfficiencyKernels::
    Chesters

Description
    Brownian coalesenceEfficiency kernel supplemented for particles in turbulent
    flows in conditions when the particle size is smaller than the local
    Kolmogorov scale.

    \f[
        P_c =
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
    Chesters.C

\*---------------------------------------------------------------------------*/

#ifndef Chesters_H
#define Chesters_H

#include "coalesenceEfficiencyKernel.H"
#include "twoPhaseSystem.H"
#include "phaseModel.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace coalesenceEfficiencyKernels
{

/*---------------------------------------------------------------------------*\
                    Class Chesters Declaration
\*---------------------------------------------------------------------------*/

class Chesters
:
    public coalesenceEfficiencyKernel
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
            const volScalarField& abscissa1,
            const volScalarField& abscissa2
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalesenceEfficiencyKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
