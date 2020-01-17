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
    Foam::populationBalanceSubModels::collisionKernels::esBGKCollision

Description
    Collision model which returns the velocity distribution to a Maxwellian
    distribution (multivariate Gaussian) with at a given collisional time
    scale. Used as a base class for esBGK as well, but a different covariance
    tensor is used as well a variable collision time scale.
    Polydispersity is also handled, but a variabel time scale is used

SourceFiles
    BGKCollision.C
    equilibriumMomentFunctions.C

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

#define defineMFunc(i,j,k)                                                      \
static void moment##i##j##k                                                     \
(                                                                               \
    mappedScalarList& Meq_,                                                     \
    const scalar& m0,                                                           \
    const scalar& u,                                                            \
    const scalar& v,                                                            \
    const scalar& w,                                                            \
    const symmTensor& sigma                                                     \
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
    typedef void (*momentFunction)
    (
        mappedScalarList& Meq,
        const scalar&,
        const scalar&,
        const scalar&,
        const scalar&,
        const symmTensor&
    );

    // Protected data

        //- Equilibrium velocity moments (for each size)
        mappedPtrList<volScalarField> Meq_;

        //- List of functions used to calculate equilibrium moments
        List<momentFunction> equilibriumMomentFunctions_;

        //- Conditioned velocity moments
        PtrList<mappedScalarList> velocityMoments_;

        //- Size conditioned equilibrium veloicity moments (Gaussian)
        PtrList<PtrList<mappedScalarList>> Gs_;

        //- Collisional times
        scalarListList Ks_;

        //- Minimum volume fraction
        scalar minM0_;


    // Protected functions

        //- Return the covariance tensor using quadrature moments
        virtual symmTensor covariance
        (
            const label celli,
            const scalar& u,
            const scalar& v,
            const scalar& w
        );

        //- Return the covariance tensor given a list of moments
        virtual symmTensor covariance
        (
            const mappedScalarList& moments,
            const scalar& u,
            const scalar& v,
            const scalar& w
        );

        //- Return the coefficient of restitution
        virtual scalar e() const
        {
            return 1.0;
        }

        //- Functions used to define equilibrium moments
        //  Stored in a list so that moment existance does not have to be 
        // checked at every cell, for every itteration. The implementation of 
        // equilibrium moment functions are in equilibriumMomentFunctions.C
        defineMFunc(0,0,0);

        defineMFunc(0,0,1);
        defineMFunc(0,1,0);
        defineMFunc(1,0,0);

        defineMFunc(0,0,2);
        defineMFunc(0,1,1);
        defineMFunc(0,2,0);
        defineMFunc(1,0,1);
        defineMFunc(1,1,0);
        defineMFunc(2,0,0);

        defineMFunc(0,0,3);
        defineMFunc(0,1,2);
        defineMFunc(0,2,1);
        defineMFunc(0,3,0);
        defineMFunc(1,0,2);
        defineMFunc(1,1,1);
        defineMFunc(1,2,0);
        defineMFunc(2,0,1);
        defineMFunc(2,1,0);
        defineMFunc(3,0,0);

        defineMFunc(0,0,4);
        defineMFunc(0,1,3);
        defineMFunc(0,2,2);
        defineMFunc(0,3,1);
        defineMFunc(0,4,0);
        defineMFunc(1,0,3);
        defineMFunc(1,3,0);
        defineMFunc(2,0,2);
        defineMFunc(2,2,0);
        defineMFunc(3,0,1);
        defineMFunc(3,1,0);
        defineMFunc(4,0,0);

        defineMFunc(0,0,5);
        defineMFunc(0,1,4);
        defineMFunc(0,2,3);
        defineMFunc(0,3,2);
        defineMFunc(0,4,1);
        defineMFunc(0,5,0);
        defineMFunc(1,0,4);
        defineMFunc(1,4,0);
        defineMFunc(2,0,3);
        defineMFunc(2,3,0);
        defineMFunc(3,0,2);
        defineMFunc(3,2,0);
        defineMFunc(4,0,1);
        defineMFunc(4,1,0);
        defineMFunc(5,0,0);

        defineMFunc(0,1,5);
        defineMFunc(0,2,4);
        defineMFunc(0,4,2);
        defineMFunc(0,5,1);
        defineMFunc(2,0,4);
        defineMFunc(2,4,0);
        defineMFunc(1,0,5);
        defineMFunc(1,5,0);
        defineMFunc(4,0,2);
        defineMFunc(4,2,0);
        defineMFunc(5,0,1);
        defineMFunc(5,1,0);

        defineMFunc(0,2,5);
        defineMFunc(0,5,2);
        defineMFunc(2,0,5);
        defineMFunc(2,5,0);
        defineMFunc(5,0,2);
        defineMFunc(5,2,0);


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

        //- Is an implicit source used
        virtual bool implicit() const
        {
            return implicit_;
        }

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

#define addMomentFunction1(map,funcs,i)                                         \
if (map.found(i))                                                               \
{                                                                               \
    funcs.append(&moment##i##00);                                               \
}

#define addMomentFunction2(map,funcs,i,j)                                       \
if (map.found(i,j))                                                             \
{                                                                               \
    funcs.append(&moment##i##j##0);                                             \
}

#define addMomentFunction3(map,funcs,i,j,k)                                     \
if (map.found(i,j,k))                                                           \
{                                                                               \
    funcs.append(&moment##i##j##k);                                             \
}

#define momentFuncHeader(i,j,k)                                                 \
void                                                                            \
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment##i##j##k\
(                                                                               \
    mappedScalarList& moments,                                                  \
    const scalar& m0,                                                           \
    const scalar& u,                                                            \
    const scalar& v,                                                            \
    const scalar& w,                                                            \
    const symmTensor& sigma                                                     \
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace collisionKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
