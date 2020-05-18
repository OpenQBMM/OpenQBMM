/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::populationBalanceSubModels::growthModels::nonLinearEvaporation

Description
    Non-linear evaporation model.

    \f[
        K_g = -\frac{E_v}{2} \left( \frac{3 V_i}{4 \pi} \right) ^{1/3}
    \f]

    where

    \vartable
        E_v         |  Evaporation coefficient [1/s]
        V_i         |  Volume of droplet i [m^3]
    \endvartable

SourceFiles
    nonLinearEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearEvaporation_H
#define nonLinearEvaporation_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                 Class nonLinearEvaporation Declaration
\*---------------------------------------------------------------------------*/

class nonLinearEvaporation
:
    public growthModel
{
public:

        //- Runtime type information
        TypeName("nonLinearEvaporation");


    // Constructors

        //- Construct from components
        nonLinearEvaporation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nonLinearEvaporation();


    // Member Functions

        //- Growth kernel
        virtual scalar Kg
        (
            const scalar& d,
            const bool lengthBased,
            const label environment
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
