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
    Foam::populationBalanceSubModels::collisionKernels::esBoltzmannCollision

Description
    Collision model which returns the velocity distribution to a Maxwellian
    distribution (multivariate Gaussian) with at a given collisional time
    scale.

SourceFiles
    BoltzmannCollision.C

\*---------------------------------------------------------------------------*/

#ifndef BoltzmannCollision_H
#define BoltzmannCollision_H

#include "collisionKernel.H"
#include "mappedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{

/*---------------------------------------------------------------------------*\
                    Class BoltzmannCollision Declaration
\*---------------------------------------------------------------------------*/

class BoltzmannCollision
:
    public collisionKernel
{
    // Private data

        //- Refrence to particle diameter
        const volScalarField& dp_;

        //- Coefficient of restitution
        scalar e_;

        //- Modified coefficient of restitution
        scalar omega_;

        //- Collisional coefficients sources
        mappedList<scalar> Is_;

        //- Collision sources
        mappedList<scalar> Cs_;


    // Protected functions

        //- Update equilibrium moments for a cell (1D)
        void I1D
        (
            const label celli,
            const label node1,
            const label node2
        );

        //- Update equilibrium moments for a cell (2D)
        void I2D
        (
            const label celli,
            const label node1,
            const label node2
        );

        //- Update equilibrium moments for a cell (3D)
        void I3D
        (
            const label celli,
            const label node1,
            const label node2
        );


public:

        //- Runtime type information
        TypeName("Boltzmann");


    // Constructors

        //- Construct from components
        BoltzmannCollision
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature,
            const bool ode
        );


    //- Destructor
    virtual ~BoltzmannCollision();


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
