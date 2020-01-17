/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
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
    Foam::firstOrderKinetic

Description
    First-order velocity moment advection with kinetic scheme.

SourceFiles
    firstOrderKineticVelocityAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef firstOrderKineticVelocityAdvection_H
#define firstOrderKineticVelocityAdvection_H

#include "velocityMomentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{

/*---------------------------------------------------------------------------*\
                    Class firstOrderKinetic Declaration
\*---------------------------------------------------------------------------*/

class firstOrderKinetic
:
    public velocityMomentAdvection
{
protected:
    // Protected data

        //- Interpolation schemes
        word weightScheme_;
        word scalarAbscissaeScheme_;
        word velocityAbscissaeScheme_;


    // Protected member functions

        // Interpolate weights and abscissae
        void interpolateNodes();



public:

    //- Runtime type information
    TypeName("firstOrderKinetic");


    // Constructors

        //- Construct from velocityMomentSet
        firstOrderKinetic
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const word& support
        );


    //- Destructor
    virtual ~firstOrderKinetic();


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
