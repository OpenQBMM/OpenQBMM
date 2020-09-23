/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::aggregationKernels::turbulentBrownian

Description
    Brownian aggregation kernel supplemented for particles in turbulent
    flows in conditions when the particle size is smaller than the local
    Kolmogorov scale.

    \f[
        K_a = \frac{2 \kappa T}{3 \mu} \frac{(\xi_i + \xi_j)^2}{\xi_i \xi_j}
            + \frac{4}{3} (\xi_i + \xi_j)^3
            \sqrt{\frac{3 \pi \varepsilon}{10 \nu}}
    \f]

    where

    \vartable
        T           |  Absolute Temperature [K]
        \varepsilon |  Turbulent dissipation rate [m2/s3]
        \kappa      |  Boltzmann's constant [J/K]
        \mu         |  Molecular dynamic iscosity [Pa s]
        \nu         |  Turbulent Kinematic viscosity [m2/s]
        \xi_i       |  Abscissa of particle i [m]
        \xi_j       |  Abscissa of particle j [m]
    \endvartable

    References
    \verbatim
        "Versuch Einer Mathematischen Theorie Der Koagulationskinetic Kolloider
        Losunger"
        M X Smoluchowski
        Zeitschrift fur Physikalische Chemie
        Volume 92, Pages 129-142, 1917
    \endverbatim

    \verbatim
        "Kinetics of turbulent coagulation studied by means of end-over-end
        rotation"
        Y Adachi, M A Cohen Stuart, R Fokkink
        Journal of Colloid and Interface Science
        Volume 165, Pages 310-317, 1994
    \endverbatim

SourceFiles
    turbulentBrownian.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentBrownian_H
#define turbulentBrownian_H

#include "aggregationKernel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class turbulentBrownian Declaration
\*---------------------------------------------------------------------------*/

class turbulentBrownian
:
    public aggregationKernel
{
    // Private data

        //- Continuous phase name
        const word continuousPhase_;

        //- Turbulence model
        const turbulenceModel& flTurb_;

        //- Temperature field
        const volScalarField& T_;

        //- Density field
        const volScalarField& rho_;

        //- Dynamic viscosity field
        const volScalarField& mu_;

        //- Turbulent dissipation rate field
        const volScalarField& epsilon_;


public:

        //- Runtime type information
        TypeName("turbulentBrownian");


    // Constructors

        //- Construct from components
        turbulentBrownian
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~turbulentBrownian();


    // Member Functions

        //- Aggregation kernel
        virtual scalar Ka
        (
            const scalar& abscissa1,
            const scalar& abscissa2,
            const vector& Ur,
            const label celli,
            const label environment = 0
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
