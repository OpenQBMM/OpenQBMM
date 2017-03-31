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
    Foam::populationBalanceSubModels::aggregationKernels::turbulentBrownian

Description
    Brownian aggregation kernel supplemented kernel for particles in turbulent
    flows in conditions when the particle size is smaller than the local
    Kolmogorov scale.

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

        //- Thermo model
        const fluidThermo& flThermo_;

        //- Turbulence model
        const compressible::turbulenceModel& flTurb_;

        //- Temperature field
        const volScalarField& T_;

        //- Kinematic viscosity field
        const volScalarField& nu_;

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
            const label& celli
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
