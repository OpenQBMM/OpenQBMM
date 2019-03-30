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

        //- Coefficient of restitution
        scalar e_;

        //- Modified coefficient of restitution
        scalar omega_;

        //- Is the Enskog term used
        Switch Enskog_;

        //- List of scalar indexes
        const labelList& scalarIndexes_;

        //- Size index
        const label sizeIndex_;

        //- List of velocity indexes
        const labelList& velocityIndexes_;

        //- Conditional velocity moment orders
        labelListList velocityMomentOrders_;

        //- Zero order collisional coefficients sources
        mappedScalarList Is_;

        //- Enskog (first order) collisional coefficients sources
        PtrList<mappedScalarList> I1s_;

        //- Collision sources
        mappedPtrList<volScalarField> Cs_;

        //- Refrence to particle diameter
        tmp<volScalarField> dp_;

        //- Refrence to particle diameter
        tmp<volScalarField> rhop_;

        //- Gradients of weights
        PtrList<volVectorField> gradWs_;

        //- Collisional fluxes
        mappedPtrList<volVectorField> Gs_;

    // Protected functions

        //- Update equilibrium moments for a cell
        void updateI
        (
            const label celli,
            const label node1,
            const label node2,
            const scalar omega
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
            const velocityQuadratureApproximation& quadrature
        );


    //- Destructor
    virtual ~BoltzmannCollision();


    // Member Functions

        //- Update unchanged fields before solving ode system
        virtual void preUpdate();

        //- Update equilibrium moments
        virtual void updateCells(const label celli);

        //- Return explicit collision source term
        virtual scalar explicitCollisionSource
        (
            const labelList& momentOrder,
            const label celli
        ) const;

        //- Return implicit collision source matrix
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volVelocityMoment& m
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace collisionKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
