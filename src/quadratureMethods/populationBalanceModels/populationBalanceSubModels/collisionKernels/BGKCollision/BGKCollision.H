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

#define defineMFunc(i,j,k) \
static void moment##i##j##k \
( \
    mappedPtrList<volScalarField>& Meq_,\
    const label,\
    const scalar&,\
    const scalar&,\
    const scalar&,\
    const scalar&,\
    const symmTensor& \
)

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

    //- Typedef of function pointer list
    typedef void (*mf)
    (
        mappedPtrList<volScalarField>& Meq_,
        const label,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const symmTensor&
    );

    // Protected data

        //- Equilibrium moments (fields)
        mappedPtrList<volScalarField> Meq_;

        //- List of functions used to calculate equilibrium moments
        List<mf> equilibriumMomentFunctions_;


    // Protected functions

        //- Return the covariance tensor
        virtual symmTensor covariance
        (
            const label celli,
            const scalar& u,
            const scalar& v,
            const scalar& w
        );

        defineMFunc(0,0,0);
        defineMFunc(1,0,0);
        defineMFunc(0,1,0);
        defineMFunc(0,0,1);
        defineMFunc(2,0,0);
        defineMFunc(0,2,0);
        defineMFunc(0,0,2);
        defineMFunc(3,0,0);
        defineMFunc(0,3,0);
        defineMFunc(0,0,3);
        defineMFunc(4,0,0);
        defineMFunc(0,4,0);
        defineMFunc(0,0,4);
        defineMFunc(5,0,0);
        defineMFunc(0,5,0);
        defineMFunc(0,0,5);
        defineMFunc(1,1,0);
        defineMFunc(1,0,1);
        defineMFunc(0,1,1);
        defineMFunc(2,1,0);
        defineMFunc(2,0,1);
        defineMFunc(0,2,1);
        defineMFunc(1,2,0);
        defineMFunc(1,0,2);
        defineMFunc(0,1,2);



public:

        //- Runtime type information
        TypeName("BGK");


    // Constructors

        //- Construct from components
        BGKCollision
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature
        );


    //- Destructor
    virtual ~BGKCollision();


    // Member Functions

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
