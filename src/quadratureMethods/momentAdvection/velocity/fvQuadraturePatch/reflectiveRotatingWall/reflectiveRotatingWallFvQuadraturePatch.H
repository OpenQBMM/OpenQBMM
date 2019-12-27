/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 Alberto Passalacqua
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
    Foam::reflectiveRotatingWallFvQuadraturePatch

Description
    Abstract class for the advection of a velocity moment set.

SourceFiles
    reflectiveRotatingWallFvQuadraturePatch.C

\*---------------------------------------------------------------------------*/

#ifndef reflectiveRotatingWallFvQuadraturePatch_H
#define reflectiveRotatingWallFvQuadraturePatch_H

#include "reflectiveFvQuadraturePatch.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class reflectiveRotatingWallFvQuadraturePatch Declaration
\*---------------------------------------------------------------------------*/

class reflectiveRotatingWallFvQuadraturePatch
:
    public reflectiveFvQuadraturePatch
{
    // Private data

        //- Origin of the rotation
        vector origin_;

        //- Axis of the rotation
        vector axis_;

        //- Rotational speed
        autoPtr<Function1<scalar>> omega_;


protected:
    // Protected data

        //- Return the wall tangential velocity
        virtual tmp<vectorField> wallTangentVelocity
        (
            const vectorField& U,
            const vectorField& n
        ) const;

public:

    //- Runtime type information
    TypeName("reflectiveRotatingWall");


    // Constructors

        //- Construct from velocityMomentSet
        reflectiveRotatingWallFvQuadraturePatch
        (
            const fvPatch& patch,
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            PtrList<surfaceVelocityNode>& nodesOwn,
            PtrList<surfaceVelocityNode>& nodesNei
        );


    //- Destructor
    virtual ~reflectiveRotatingWallFvQuadraturePatch();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
