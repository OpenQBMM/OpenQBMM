/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::noAdvection

Description
    Empty scheme for cases without moment advection.

SourceFiles
    noAdvection.C
    noAdvectionI.H

\*---------------------------------------------------------------------------*/

#ifndef noAdvection_H
#define noAdvection_H

#include "velocityMomentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{

/*---------------------------------------------------------------------------*\
                        Class noAdvection Declaration
\*---------------------------------------------------------------------------*/

class noAdvection
:
    public velocityMomentAdvection
{
public:

    //- Runtime type information
    TypeName("noAdvection");


    // Constructors

        //- Construct from velocityMomentSet
        noAdvection
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const word& support
        );


    //- Destructor
    virtual ~noAdvection();


    // Public member functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo();

        //- Update moment advection
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace velocityAdvection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
