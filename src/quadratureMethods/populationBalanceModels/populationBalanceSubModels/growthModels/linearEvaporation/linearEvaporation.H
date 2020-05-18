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
    Foam::populationBalanceSubModels::growthModels::linearEvaporation

Description
    Linear evaporation model.

    \f[
        C_g = -E_v \frac{\pi d_i^3}{6}
    \f]

    where

    \vartable
        E_v         |  Evaporation coefficient [m^2/s]
        d_i         |  Diameter of droplet i [m^3]
    \endvartable

SourceFiles
    linearEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef linearEvaporation_H
#define linearEvaporation_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class linear Declaration
\*---------------------------------------------------------------------------*/

class linearEvaporation
:
    public growthModel
{
public:

        //- Runtime type information
        TypeName("linearEvaporation");


    // Constructors

        //- Construct from components
        linearEvaporation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~linearEvaporation();


    // Member Functions

        //- Aggregation kernel
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
