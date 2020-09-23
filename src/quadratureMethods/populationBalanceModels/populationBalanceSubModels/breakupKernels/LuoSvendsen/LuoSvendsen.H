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
    Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen

Description
    Breakup kernel for liquid-solid systems.

    Reference
    \verbatim
        "Theoretical model for drop and bubble breakup in turbulent dispersion"
        H Luo, H F Svendsen
        AIChE Journal
        Volume 42, Pages 1225-1233, 1996
    \endverbatim

SourceFiles
    LuoSvendsen.C

\*---------------------------------------------------------------------------*/

#ifndef LuoSvendsen_H
#define LuoSvendsen_H

#include "breakupKernel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{

/*---------------------------------------------------------------------------*\
                    Class LuoSvendsen Declaration
\*---------------------------------------------------------------------------*/

class LuoSvendsen
:
    public breakupKernel
{
private:

    // Private data

        //- Continuous phase name
        const word continuousPhase_;

        //- Model coefficient
        dimensionedScalar Cb_;

        //- Exponent of the turbulent dissipation rate
        scalar epsilonExp_;

        //- Exponent of the kinematic viscosity
        scalar nuExp_;

        //- Exponent of the abscissa
        scalar sizeExp_;

        //- Turbulence model
        const turbulenceModel& flTurb_;

        //- Turbulent dissipation rate
        const volScalarField& epsilon_;

        //- Dynamic viscosity
        const volScalarField& mu_;

        //- Density
        const volScalarField& rho_;

public:

        //- Runtime type information
        TypeName("LuoSvendsen");


    // Constructors

        //- Construct from components
        LuoSvendsen
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~LuoSvendsen();


    // Member Functions

        //- Breakup kernel
        virtual scalar Kb
        (
            const scalar& abscissa,
            const label celli,
            const label environment = 0
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
