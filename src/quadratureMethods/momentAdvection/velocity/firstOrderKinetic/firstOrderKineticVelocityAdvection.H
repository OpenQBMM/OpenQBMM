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
    Foam::firstOrderKineticVelocityAdvection

Description
    First-order velocity moment advection with kinetic scheme.

SourceFiles
    firstOrderKinetic.C
    firstOrderKineticI.H

\*---------------------------------------------------------------------------*/

#ifndef firstOrderKinetic_H
#define firstOrderKinetic_H

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
    // Private data

        //- List of quadrature nodes used for moment advection
        const mappedPtrList<volVectorNode>& nodes_;

        //- List of interpolated nodes (neighbour)
        autoPtr<mappedPtrList<surfaceVectorNode>> nodesNei_;

        //- List of interpolated nodes (owner)
        autoPtr<mappedPtrList<surfaceVectorNode>> nodesOwn_;

        //- Flux moment inversion method
//         autoPtr<hyperbolicFieldMomentInversion> momentFieldInverter_;


protected:

        //- Protected member functions
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
