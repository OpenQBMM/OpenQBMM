/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::aggregationKernels::shearRate

Description
    Coalescence model for droplets of O'Rourke

    \f[
        \omega = \pi \left( d_i + d_j \right)^2 |U_i - U_j|
    \f]

    where

    \vartable
        d_i         |  Diameter of bubble i [m]
        d_j         |  Diameter of bubble j [m]
        U_i         |  Velocity of bubble i [m]
        U_j         |  Velocity of bubble j [m]
    \endvartable

    References
    \verbatim

    \endverbatim

SourceFiles
    shearRateFrequency.C

\*---------------------------------------------------------------------------*/

#ifndef shearRateFrequency_H
#define shearRateFrequency_H

#include "coalescenceFrequencyKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace coalescenceFrequencyKernels
{

/*---------------------------------------------------------------------------*\
                    Class shearRate Declaration
\*---------------------------------------------------------------------------*/

class shearRate
:
    public coalescenceFrequencyKernel
{
public:

        //- Runtime type information
        TypeName("shearRate");


    // Constructors

        //- Construct from components
        shearRate
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& continuousPhase
        );


    //- Destructor
    virtual ~shearRate();


    // Member Functions

        //- Aggregation kernel
        virtual scalar omega
        (
            const scalar& abscissa1,
            const scalar& abscissa2,
            const vector& Ur,
            const label celli
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceFrequencyKernels
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
