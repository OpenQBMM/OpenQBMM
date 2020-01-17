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
    Copyright (C) 2019 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::breakupKernels::exponentialBreakup

Description
    Exponential breakup kernel. The value of the kernel is read from dictionary
    or assumed equal to 1. The minimum value of the abscissa after which the
    model is applied can be defined by the user or it is assumed to be 1.

SourceFiles
    exponentialBreakup.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialBreakup_H
#define exponentialBreakup_H

#include "breakupKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{

/*---------------------------------------------------------------------------*\
                    Class exponentialBreakup Declaration
\*---------------------------------------------------------------------------*/

class exponentialBreakup
:
    public breakupKernel
{
private:

    // Private data

        //- Minimum value of the abscissa to apply the model
        dimensionedScalar minAbscissa_;

        //- Exponent of the kernel
        dimensionedScalar expCoeff_;

public:

        //- Runtime type information
        TypeName("exponential");


    // Constructors

        //- Construct from components
        exponentialBreakup
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~exponentialBreakup();


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
