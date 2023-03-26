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
    Foam::populationBalanceSubModels::growthModels::LageCondensation

Description
    Condensation model of Lage (2011). 

    Reference
    \verbatim
        "On the representation of QMOM as aweighted-residual method — 
        The dual-quadrature method of generalized moments"
        P L C Lage
        Computers & Chemical Engineering
        Volume 35, Issue 11, Pages 2186-2203, 2011
    \endverbatim

SourceFiles
    LageCondensation.C

\*---------------------------------------------------------------------------*/

#ifndef LageCondensation_H
#define LageCondensation_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class LageCondensation Declaration
\*---------------------------------------------------------------------------*/

class LageCondensation
:
    public growthModel
{
public:

        //- Runtime type information
        TypeName("LageCondensation");


    // Constructors

        //- Construct from components
        LageCondensation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~LageCondensation();


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
