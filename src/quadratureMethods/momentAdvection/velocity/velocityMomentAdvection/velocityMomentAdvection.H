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
#include "quadratureApproximations.H"
#include "surfaceInterpolationScheme.H"
#include "fvQuadraturePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class velocityMomentAdvection Declaration
\*---------------------------------------------------------------------------*/

class velocityMomentAdvection
{
protected:

    // Protected data

        //- Name of the distribution associated to the quadrature approximation
        const word name_;

        //- Reference to the moments to advect
        const volVelocityMomentFieldSet& moments_;

        //- List of quadrature nodes used for moment advection
        const PtrList<volVelocityNode>& nodes_;

        //- List of interpolated nodes (neighbour)
        autoPtr<PtrList<surfaceVelocityNode>> nodesNei_;

        //- List of interpolated nodes (owner)
        autoPtr<PtrList<surfaceVelocityNode>> nodesOwn_;

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

        //- List of boundaries conditions
        PtrList<fvQuadraturePatch> boundaries_;


    //- Protected member function

            //- Update fluxes due to wall callisions
            void updateBoundaryConditions();


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

        //- Disallow default bitwise copy construct
        velocityMomentAdvection(const velocityMomentAdvection&) = delete;


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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const velocityMomentAdvection&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "velocityMomentAdvectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
