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
    Foam::reflectiveFvQuadraturePatch

Description
    Abstract class for the advection of a velocity moment set.

SourceFiles
    reflectiveFvQuadraturePatch.C

\*---------------------------------------------------------------------------*/

#ifndef reflectiveFvQuadraturePatch_H
#define reflectiveFvQuadraturePatch_H

#include "fvQuadraturePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class reflectiveFvQuadraturePatch Declaration
\*---------------------------------------------------------------------------*/

class reflectiveFvQuadraturePatch
:
    public fvQuadraturePatch
{
    // Private data

        //- Coefficient of restitution
        scalar ew_;

public:

    //- Runtime type information
    TypeName("reflective");


    // Constructors

        //- Construct from velocityMomentSet
        reflectiveFvQuadraturePatch
        (
            const fvPatch& patch,
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            PtrList<surfaceVelocityNode>& nodesOwn,
            PtrList<surfaceVelocityNode>& nodesNei
        );


    //- Destructor
    virtual ~reflectiveFvQuadraturePatch();

    // Member Functions

        //- Update boundary flux
        virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //