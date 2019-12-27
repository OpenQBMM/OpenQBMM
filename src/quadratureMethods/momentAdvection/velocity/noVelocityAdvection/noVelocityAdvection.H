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
    Foam::noAdvection

Description
    Empty scheme for cases without moment advection. Only to be used for
    zero-dimensional cases.

SourceFiles
    noAdvectionVelocityAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef noAdvectionVelocityAdvection_H
#define noAdvectionVelocityAdvection_H

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
        virtual scalar realizableCo() const;

        //- Return the maximum Courant number based on velocities
        virtual scalar CoNum() const;

        //- Update moment advection
        virtual void update();

        //- Update moment advection using a specified velocity field
        virtual void update
        (
            const surfaceScalarField& U,
            const bool wallCollisions
        );

        //- Update moment advection using a specified list of velocity fields
        virtual void update
        (
            const mappedPtrList<volVectorField>& Us,
            const bool wallCollisions
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace velocityAdvection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
