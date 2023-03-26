/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::growthModels::PruppacherKlett

Description
    Diffusion-driven growth model. The typical value for Cg is 0.78.
    
    Reference
    \verbatim
        "Microphysics of Clouds and Precipitation"
        H R Pruppacher, J D Klett
        Reidel, Boston, MA, 1980
    \endverbatim

SourceFiles
    diffusionControlledGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionControlledGrowth_H
#define diffusionControlledGrowth_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
            Class diffusionControlledGrowth Declaration
\*---------------------------------------------------------------------------*/

class diffusionControlledGrowth
:
    public growthModel
{
public:

        //- Runtime type information
        TypeName("diffusionControlledGrowth");


    // Constructors

        //- Construct from components
        diffusionControlledGrowth
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~diffusionControlledGrowth();


    // Member Functions

        //- Growth model
        virtual scalar Kg
        (
            const scalar& abscissa,
            const bool lengthBased,
            const label environment = 0
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
