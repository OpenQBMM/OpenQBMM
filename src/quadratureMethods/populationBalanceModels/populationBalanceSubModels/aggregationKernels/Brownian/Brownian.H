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
    Foam::populationBalanceSubModels::aggregationKernels::Brownian

Description
    Brownian aggregation kernel.

    \f[
        K_a = \frac{2 \kappa T}{3 \mu} \frac{(\xi_i + \xi_j)^2}{\xi_i \xi_j}
    \f]

    where

    \vartable
        T           |  Absolute Temperature [K]
        \kappa      |  Boltzmann's constant [J/K]
        \mu         |  Molecular Viscosity [Pa s]
        \xi_i       |  Abscissa of particle i [m]
        \xi_j       |  Abscissa of particle j [m]
    \endvartable

    Reference
    \verbatim
        "Versuch Einer Mathematischen Theorie Der Koagulationskinetic Kolloider
        Losunger"
        M X Smoluchowski
        Zeitschrift fur Physikalische Chemie
        Volume 92, Pages 129-142, 1917
    \endverbatim

SourceFiles
    Brownian.C

\*---------------------------------------------------------------------------*/

#ifndef Brownian_H
#define Brownian_H

#include "aggregationKernel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class Brownian Declaration
\*---------------------------------------------------------------------------*/

class Brownian
:
    public aggregationKernel
{
    //- Private data

        //- Continuous phase name
        const word continuousPhase_;

        //- Temperature field
        const volScalarField& T_;

        //- Viscosity field
        const volScalarField& mu_;


public:

        //- Runtime type information
        TypeName("Brownian");


    // Constructors

        //- Construct from components
        Brownian
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Brownian();


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
