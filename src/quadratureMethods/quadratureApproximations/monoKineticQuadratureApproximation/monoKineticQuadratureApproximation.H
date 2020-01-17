/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
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
    Foam::monoKineticQuadratureApproximation

Description
    Uses quadrature approximation as a base to handle size moments, while the
    velocity moments and abscissae are stored in this class. The one
    dimensional CQMOM inversion is also done.

SourceFiles
    monoKineticQuadratureApproximation.C
    monoKineticQuadratureApproximationI.H

\*---------------------------------------------------------------------------*/

#ifndef monoKineticQuadratureApproximation_H
#define monoKineticQuadratureApproximation_H

#include "fvCFD.H"
#include "quadratureApproximation.H"
#include "PtrList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "univariateMomentSet.H"
#include "quadratureNodes.H"
#include "momentFieldSets.H"
#include "fieldMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class monoKineticQuadratureApproximation Declaration
\*---------------------------------------------------------------------------*/

class monoKineticQuadratureApproximation
:
    public quadratureApproximation<volScalarMoment, volScalarNode>
{
    // Private data

        //- Refrence to mean velocity
        const volVectorField& U_;

        //- Number of primary quadrature nodes
        label nNodes_;

        //- List of velocity moments
        PtrList<volVectorField> velocityMoments_;

        //- List of velocity abscissae
        PtrList<volVectorField> velocityAbscissae_;

        //- List of interpolated nodes (neighbour)
        autoPtr<PtrList<surfaceScalarNode>> nodesNei_;

        //- List of interpolate velocities (neighbour)
        PtrList<surfaceVectorField> velocitiesNei_;

        //- List of interpolated nodes (owner)
        autoPtr<PtrList<surfaceScalarNode>> nodesOwn_;

        //- List of interpolate velocities (owner)
        PtrList<surfaceVectorField> velocitiesOwn_;

        //- Minimum values for m1 and m0
        //  minM1 should be the same order as residual alpha
        scalar minM0_;
        scalar minM1_;


public:

    // Constructors

        //- Construct from mesh
        explicit monoKineticQuadratureApproximation
        (
            const word& name,
            const fvMesh& mesh,
            const word& support
        );


    //- Destructor
    virtual ~monoKineticQuadratureApproximation();


    // Member Functions

        // Access

            //- Const access to the velocity abscissae
            inline const PtrList<volVectorField>& velocities() const;

            //- Const access to the velocity abscissae
            inline PtrList<volVectorField>& velocities();

            //- Const access to the velocity moments
            inline const PtrList<volVectorField>& velocityMoments() const;

            //- Non-const access to the velocity moments
            inline PtrList<volVectorField>& velocityMoments();

            //- Const access to the interpolated moments (owner)
            inline const PtrList<surfaceScalarNode>& nodesOwn() const;

            //- Const access to the interpolated velocities (owner)
            inline const PtrList<surfaceVectorField>& velocitiesOwn() const;

            //- Const access to the interpolated moments (neighbour)
            inline const PtrList<surfaceScalarNode>& nodesNei() const;

            //- Const access to the interpolated velocities (neighbour)
            inline const PtrList<surfaceVectorField>& velocitiesNei() const;


        // Edit

            //- Interpolate nodes
            void interpolateNodes();

            //- Update node velocities on fixed-value boundaries
            void updateBoundaryVelocities();

            //- Recalculate the velocity abscissae from the velocity moments
            void updateVelocities();

            //- Recalculate the velocity abscissae from the velocity moments
            void updateLocalVelocities(const label celli);

            //- Recalculate the quadrature nodes from the moments
            void updateAllQuadrature();

            //- Recalculate the quadrature nodes from the moments
            bool updateAllLocalQuadrature
            (
                const label celli,
                const bool failOnRealizability = true
            );

            //- Recalculate the moments from the quadrature nodes
            void updateVelocityMoments();

            //- Recalculate the moments from the quadrature nodes
            void updateLocalVelocityMoments(const label celli);

            //- Recalculate all moments from the quadrature nodes
            void updateAllMoments();

            //- Recalculate all moments from the quadrature nodes
            void updateAllLocalMoments(const label celli);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "monoKineticQuadratureApproximationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
