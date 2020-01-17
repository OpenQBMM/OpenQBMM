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
    Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision

Description
    Analytical expansion of Boltzmann collision intergral for monodisperse and
    polydisperse particulate systems

    Reference:
    \verbatim
        "Computational Models for Polydisperse Particulate and Multiphase Systems"
        Marchisio, D.L., Fox, R.O.,
        Cambridge Core, 2013
    \endverbatim

SourceFiles
    BoltzmannCollision.C
    ICoefficientFunctions.C
    IxCoefficientFunctions.C
    IyCoefficientFunctions.C
    IzCoefficientFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef BoltzmannCollision_H
#define BoltzmannCollision_H

#include "collisionKernel.H"
#include "mappedLists.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{

#define defineIFuncs(i,j,k)                                                     \
static void I##i##j##k                                                          \
(                                                                               \
    mappedScalarList& Is,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
);                                                                              \
                                                                                \
static void Ix##i##j##k                                                         \
(                                                                               \
    mappedScalarList& Is,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
);                                                                              \
                                                                                \
static void Iy##i##j##k                                                         \
(                                                                               \
    mappedScalarList& Is,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
);                                                                              \
                                                                                \
static void Iz##i##j##k                                                         \
(                                                                               \
    mappedScalarList& Is,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
);

/*---------------------------------------------------------------------------*\
                    Class BoltzmannCollision Declaration
\*---------------------------------------------------------------------------*/

class BoltzmannCollision
:
    public collisionKernel
{
    //- Typedef of function pointer list
    typedef void (*momentFunction)
    (
        mappedScalarList& Is,
        const scalarList& omegaPow,
        const vectorList& gPow,
        const scalar& gMagSqr,
        const vectorList& vPow
    );

    // Private data

        //- Coefficient of restitution
        scalar e_;

        //- Modified coefficient of restitution
        scalar omega_;

        //- Is the Enskog term used
        Switch Enskog_;

        //- List of scalar indexes
        const labelList& scalarIndexes_;

        //- Zero order collisional coefficients sources
        mappedScalarList Is_;

        //- List of functions used to calculate equilibrium moments
        List<momentFunction> coefficientFunctions_;

        //- Enskog (first order) collisional coefficients sources
        PtrList<mappedScalarList> I1s_;

        //- List of functions used to calculate equilibrium moments
        List<List<momentFunction>> enskogFunctions_;

        //- Collision sources
        mappedPtrList<volScalarField> Cs_;

        //- Gradients of weights
        PtrList<volVectorField> gradWs_;

        //- Collisional fluxes
        mappedPtrList<volVectorField> Gs_;

    // Protected functions

        //- Functions used to define integral coefficients
        defineIFuncs(0,0,0)

        defineIFuncs(0,0,1)
        defineIFuncs(0,1,0)
        defineIFuncs(1,0,0)

        defineIFuncs(0,0,2)
        defineIFuncs(0,1,1)
        defineIFuncs(0,2,0)
        defineIFuncs(1,0,1)
        defineIFuncs(1,1,0)
        defineIFuncs(2,0,0)

        defineIFuncs(0,0,3)
        defineIFuncs(0,1,2)
        defineIFuncs(0,2,1)
        defineIFuncs(0,3,0)
        defineIFuncs(1,0,2)
        defineIFuncs(1,1,1)
        defineIFuncs(1,2,0)
        defineIFuncs(2,0,1)
        defineIFuncs(2,1,0)
        defineIFuncs(3,0,0)

        defineIFuncs(0,0,4)
//         defineIFuncs(0,1,3)
//         defineIFuncs(0,2,2)
//         defineIFuncs(0,3,1)
        defineIFuncs(0,4,0)
//         defineIFuncs(1,0,3)
//         defineIFuncs(1,3,0)
//         defineIFuncs(2,0,2)
//         defineIFuncs(2,2,0)
//         defineIFuncs(3,0,1)
//         defineIFuncs(3,1,0)
        defineIFuncs(4,0,0)


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

#define addIFunction1(map,i)                                                    \
if (map.found(i))                                                               \
{                                                                               \
    coefficientFunctions_.append(&I##i##00);                                    \
    if (Enskog_)                                                                \
    {                                                                           \
        enskogFunctions_[0].append(&Ix##i##00);                                 \
        if (velocityIndexes_.size() > 1)                                        \
        {                                                                       \
            enskogFunctions_[1].append(&Iy##i##00);                             \
        }                                                                       \
        if (velocityIndexes_.size() > 2)                                        \
        {                                                                       \
            enskogFunctions_[2].append(&Iz##i##00);                             \
        }                                                                       \
    }                                                                           \
}

#define addIFunction2(map,i,j)                                                  \
if (map.found(i,j))                                                             \
{                                                                               \
    coefficientFunctions_.append(&I##i##j##0);                                  \
    if (Enskog_)                                                                \
    {                                                                           \
        enskogFunctions_[0].append(&Ix##i##j##0);                               \
        if (velocityIndexes_.size() > 1)                                        \
        {                                                                       \
            enskogFunctions_[1].append(&Iy##i##j##0);                           \
        }                                                                       \
        if (velocityIndexes_.size() > 2)                                        \
        {                                                                       \
            enskogFunctions_[2].append(&Iz##i##j##0);                           \
        }                                                                       \
    }                                                                           \
}

#define addIFunction3(map,i,j,k)                                                \
if (map.found(i,j,k))                                                           \
{                                                                               \
    coefficientFunctions_.append(&I##i##j##k);                                  \
    if (Enskog_)                                                                \
    {                                                                           \
        enskogFunctions_[0].append(&Ix##i##j##k);                               \
        if (velocityIndexes_.size() > 1)                                        \
        {                                                                       \
            enskogFunctions_[1].append(&Iy##i##j##k);                           \
        }                                                                       \
        if (velocityIndexes_.size() > 2)                                        \
        {                                                                       \
            enskogFunctions_[2].append(&Iz##i##j##k);                           \
        }                                                                       \
    }                                                                           \
}

#define IFuncHeader(i,j,k)                                                      \
void                                                                            \
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::        \
I##i##j##k                                                                      \
(                                                                               \
    mappedScalarList& Is,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
)

#define IxFuncHeader(i,j,k)                                                     \
void                                                                            \
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::        \
Ix##i##j##k                                                                     \
(                                                                               \
    mappedScalarList& Ix,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
)

#define IyFuncHeader(i,j,k)                                                     \
void                                                                            \
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::        \
Iy##i##j##k                                                                     \
(                                                                               \
    mappedScalarList& Iy,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
)

#define IzFuncHeader(i,j,k)                                                     \
void                                                                            \
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::        \
Iz##i##j##k                                                                     \
(                                                                               \
    mappedScalarList& Iz,                                                       \
    const scalarList& omegaPow,                                                 \
    const vectorList& gPow,                                                     \
    const scalar& gMagSqr,                                                  \
    const vectorList& vPow                                                      \
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace collisionKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
