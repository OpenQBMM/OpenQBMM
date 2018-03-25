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
    Foam::velocityMomentAdvection

Description
    Abstract class for the advection of a velocity moment set.

SourceFiles
    velocityMomentAdvection.C
    velocityMomentAdvectionI.H

\*---------------------------------------------------------------------------*/

#ifndef velocityMomentAdvection_H
#define velocityMomentAdvection_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "quadratureNodes.H"
#include "momentFieldSets.H"
#include "hyperbolicFieldMomentInversion.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class velocityMomentAdvection Declaration
\*---------------------------------------------------------------------------*/

class velocityMomentAdvection
{
    // Private member functions

        //- Disallow default bitwise copy construct
        velocityMomentAdvection(const velocityMomentAdvection&);

        //- Disallow default bitwise assignment
        void operator=(const velocityMomentAdvection&);


protected:

    // Protected data

        //- Name of the distribution associated to the quadrature approximation.
        const word name_;

        //- Reference to the moments to advect
        const volVectorMomentFieldSet& moments_;

        //- Number of moments
        const label nMoments_;

        //- Field used to find upwind values of the owner
        surfaceScalarField own_;

        //- Field used to find upwind values of the neighbour
        surfaceScalarField nei_;

        //- Support of the distribution function
        word support_;

        //- Moment indicies ues to map moments
        labelListList momentOrders_;

        //- Node indicies used to map nodes
        labelListList nodeIndexes_;

        //- Advection terms for moment transport equations
        PtrList<volScalarField> divMoments_;


public:

    //- Runtime type information
    TypeName("velocityMomentAdvection");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            velocityMomentAdvection,
            dictionary,
            (
                const dictionary& dict,
                const velocityQuadratureApproximation& quadrature,
                const word& support
            ),
            (dict, quadrature, support)
        );


    // Constructors

        //- Construct from velocityMomentSet
        velocityMomentAdvection
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const word& support
        );


    //- Destructor
    virtual ~velocityMomentAdvection();


    // Selectors

        static autoPtr<velocityMomentAdvection> New
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const word& support
        );


    // Member Functions

        //- Return the advection term for the moment transport equation
        inline const PtrList<volScalarField>& divMoments() const;

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const = 0;

        //- Return the maximum Courant number based on velocities
        virtual scalar CoNum() const = 0;

        //- Update moment advection
        virtual void update() = 0;

        //- Update moment advection using a specified velocity field
        virtual void update
        (
            const surfaceScalarField& U,
            const bool wallCollisions
        ) = 0;

        //- Update moment advection using a specified list of velocity fields
        virtual void update
        (
            const mappedPtrList<volVectorField>& Us,
            const bool wallCollisions
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "velocityMomentAdvectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
