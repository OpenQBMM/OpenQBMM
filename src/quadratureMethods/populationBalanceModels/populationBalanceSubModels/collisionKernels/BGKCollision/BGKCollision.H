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
    Foam::populationBalanceSubModels::collisionKernels::esBGKCollision

Description
    Collision model which returns the velocity distribution to a Maxwellian
    distribution (multivariate Gaussian) with at a given collisional time
    scale.

SourceFiles
    BGKCollision.C

\*---------------------------------------------------------------------------*/

#ifndef BGKCollision_H
#define BGKCollision_H

#include "collisionKernel.H"
#include "mappedPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{

/*---------------------------------------------------------------------------*\
                    Class BGKCollision Declaration
\*---------------------------------------------------------------------------*/

class BGKCollision
:
    public collisionKernel
{
    // Private data

        //- Collisional time scale
        dimensionedScalar tauCollisional_;


protected:
    // Protected data

        //- Equilibrium moments (fields)
        mappedPtrList<volScalarField> Meqf_;

        //- Equilibrium moments (fields)
        mappedList<scalar> Meq_;


    // Protected functions

        //- Update equilibrium moments for a cell (1D)
        void updateCells1D(const label celli);

        //- Update equilibrium moments for a cell (2D)
        void updateCells2D(const label celli);

        //- Update equilibrium moments for a cell (3D)
        void updateCells3D(const label celli);

        //- Update equilibrium moments for fields (1D)
        void updateFields1D();

        //- Update equilibrium moments for fields (2D)
        void updateFields2D();

        //- Update equilibrium moments for fields (3D)
        void updateFields3D();


public:

        //- Runtime type information
        TypeName("BGK");


    // Constructors

        //- Construct from components
        BGKCollision
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature,
            const bool ode
        );


    //- Destructor
    virtual ~BGKCollision();


    // Member Functions

        //- Update equilibrium moments
        virtual void updateCells(const label celli);

        //- Update equilibrium moments
        virtual void updateFields();

        //- Return explicit collision source term
        virtual scalar explicitCollisionSource
        (
            const label mi,
            const label celli
        ) const;

        //- Return implicit collision source matrix
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volVectorMoment& m
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace collisionKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
