/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::populationBalanceSubModels::breakupKernels::tubBreakup

Description
    Turbulent breakup kernel.

SourceFiles
    tubBreakup.C

\*---------------------------------------------------------------------------*/

#ifndef tubBreakup_H
#define tubBreakup_H

#include "breakupKernel.H"
#include "volFields.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class turbulenceModel;

namespace populationBalanceSubModels
{
namespace breakupKernels
{

/*---------------------------------------------------------------------------*\
                           Class tubBreakup Declaration
\*---------------------------------------------------------------------------*/

class tubBreakup
:
    public breakupKernel
{
    // Private Data

        //- Continuous phase name
        word continuousPhase_;

        //- Density
        dimensionedScalar des_;

        //- Viscosity
        dimensionedScalar vis_;

        //- Fractal dimension
        dimensionedScalar Df_;

        //- Initial radius
        dimensionedScalar R0_;

        //- Energy dissipation rate
        //dimensionedScalar eps_;

        //- Flag to use constant epsilon
        bool useConstantEpsilon_;

        //- Constant epsilon value
        dimensionedScalar constantEpsilon_;

        //- Pointer to the turbulence model
        const turbulenceModel* flTurb_;

        //- Epsilon field
        tmp<volScalarField> epsilon_;
volScalarField gField_;

public:

    //- Runtime type information
    TypeName("tubBreakup");


    // Constructors

        //- Construct from components
        tubBreakup
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~tubBreakup();


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
