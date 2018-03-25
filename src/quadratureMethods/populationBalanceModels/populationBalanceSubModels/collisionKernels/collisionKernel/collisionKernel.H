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
    Foam::populationBalanceSubModels::collisionKernel

Description
    Abstract class for collision kernels.

SourceFiles
    collisionKernel.C

\*---------------------------------------------------------------------------*/

#ifndef collisionKernel_H
#define collisionKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "quadratureApproximations.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                    Class collisionKernel Declaration
\*---------------------------------------------------------------------------*/

class collisionKernel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        collisionKernel(const collisionKernel&);

        //- Disallow default bitwise assignment
        void operator=(const collisionKernel&);


protected:

    // Protected data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Quadrature approximation
        const velocityQuadratureApproximation& quadrature_;

        //- Moment indicies ues to map moments
        const labelListList momentOrders_;

        //- Number of dimensions
        const label nDimensions_;


    // Protected functions

        //- Variances for a cell
        scalar meanVelocity(const scalar& m0, const scalar& m1) const;

        //- Mean velocities for a field
        tmp<volScalarField> meanVelocity
        (
            const volScalarField& m0,
            const volScalarField& m1
        );

        //- Variances for a cell
        scalar variance
        (
            const scalar& m0,
            const scalar& sqrM1,
            const scalar& m2
        ) const;

        //- Variance for a field
        tmp<volScalarField> variance
        (
            const volScalarField& m0,
            const volScalarField& sqrM1,
            const volScalarField& m2
        );


public:

    //- Runtime type information
    TypeName("collisionKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        collisionKernel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature,
            const bool ode
        ),
        (dict, mesh, quadrature, ode)
    );


    // Constructors

        //- Construct from components
        collisionKernel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature,
            const bool ode
        );


    // Selectors

        static autoPtr<collisionKernel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature,
            const bool ode
        );


    //- Destructor
    virtual ~collisionKernel();


    // Member Functions

        //- Update collision sources in a cell
        virtual void updateCells(const label celli) = 0;

        //- Update collision sources
        virtual void updateFields() = 0;

        //- Return explicit collision source term
        virtual scalar explicitCollisionSource
        (
            const label mi,
            const label celli
        ) const = 0;

        //- Return implicit collision source matrix
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volVectorMoment& m
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
