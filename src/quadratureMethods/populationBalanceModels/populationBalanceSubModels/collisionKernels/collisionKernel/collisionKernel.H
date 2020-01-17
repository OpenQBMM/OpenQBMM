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
protected:

    // Protected data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Quadrature approximation
        const velocityQuadratureApproximation& quadrature_;

        //- Moment indicies ues to map moments
        const labelListList& momentOrders_;

        //- Nodes indexes
        const labelListList& nodeIndexes_;

        //- Velocity indexes
        const labelList& velocityIndexes_;

        //- Number of dimensions
        const label nDimensions_;

        //- List of velocity moment orders
        labelListList velocityMomentOrders_;

        //- List of velocity node indexes
        labelListList velocityNodeIndexes_;

        //- Size index
        label sizeIndex_;

        //- Number of sizes
        label nSizes_;

        //- Refrence to particle diameter
        tmp<volScalarField> dp_;

        //- List of patricle densities
        scalarList rhos_;

        //- Minimum diameter
        scalar minD_;

        //- Is an implicit solution used
        Switch implicit_;


    // Protected functions

        static tmp<volScalarField> lookupOrInitialize
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& dict,
            const word& entryName,
            const dimensionSet& dims
        );

        //- Return the diameter given a node and cell
        scalar d(const label node, const label celli) const;

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
            const velocityQuadratureApproximation& quadrature
        ),
        (dict, mesh, quadrature)
    );


    // Constructors

        //- Construct from components
        collisionKernel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature
        );

        //- Disallow default bitwise copy construct
        collisionKernel(const collisionKernel&) = delete;


    // Selectors

        static autoPtr<collisionKernel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature
        );


    //- Destructor
    virtual ~collisionKernel();


    // Member Functions

        //- Is an implicit source used
        virtual bool implicit() const
        {
            return implicit_;
        }

        //- Update unchanged fields before solving ode system
        virtual void preUpdate();

        //- Update collision sources in a cell
        virtual void updateCells(const label celli) = 0;

        //- Update collision sources
        virtual void updateFields();

        //- Return explicit collision source term
        virtual scalar explicitCollisionSource
        (
            const labelList& momentOrder,
            const label celli
        ) const = 0;

        //- Return implicit collision source matrix
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volVelocityMoment& m
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const collisionKernel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
