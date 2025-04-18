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
    Copyright (C) 2019-2023 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::aggregationKernels::differentialForceAggregation

Description
    Differential force aggregation kernel with user defined coefficient
    (default Cagg = 1).

    \f[
        K_a = C_a (\xi_i + \xi_j)^2 |\xi_i^2 - \xi_j^2|
    \f]

    where

    \vartable
        C_a         |  Aggregation coefficient [1/s]
        \xi_i       |  Abscissa of particle i [m]
        \xi_j       |  Abscissa of particle j [m]
    \endvartable

SourceFiles
    differentialForceAggregation.C

\*---------------------------------------------------------------------------*/

#ifndef differentialForceAggregation_H
#define differentialForceAggregation_H

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                 Class differentialForceAggregation Declaration
\*---------------------------------------------------------------------------*/

class differentialForceAggregation
:
    public aggregationKernel
{
public:

        //- Runtime type information
        TypeName("differentialForce");


    // Constructors

        //- Construct from components
        differentialForceAggregation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~differentialForceAggregation();


    // Member Functions

        //- Aggregation kernel
        virtual scalar Ka
        (
            const scalar& d1,
            const scalar& d2,
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
