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

\*---------------------------------------------------------------------------*/

#include "BoltzmannCollision.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{
    defineTypeNameAndDebug(BoltzmannCollision, 0);

    addToRunTimeSelectionTable
    (
        collisionKernel,
        BoltzmannCollision,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

void
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::updateI
(
    const label celli,
    const label node1,
    const label node2,
    const scalar omega
)
{
    const vector& u1 = quadrature_.nodes()[node1].velocityAbscissae()[celli];
    const vector& u2 = quadrature_.nodes()[node2].velocityAbscissae()[celli];

    // Store powers of velocities to reduce evaluations
    vector g = u1 - u2;
    scalarList omegaPow(6, omega);
    vectorList gPow(6, g);
    scalar gMagSqr(magSqr(g));
    vectorList vPow(6, u1);

    forAll(gPow, powi)
    {
        omegaPow[powi] = pow(omega, powi);
        for (label cmpt = 0; cmpt < 3; cmpt++)
        {
            gPow[powi][cmpt] = pow(g[cmpt], powi);
            vPow[powi][cmpt] = pow(u1[cmpt], powi);
        }
    }

    // Update coefficients for zero order terms
    forAllIter(List<momentFunction>, coefficientFunctions_, iter)
    {
        (*iter)(Is_, omegaPow, gPow, gMagSqr, vPow);
    }

    // Update first order (Enskog) terms if used
    if (Enskog_)
    {
        forAllIter(List<momentFunction>, enskogFunctions_[0], iter)
        {
            (*iter)(I1s_[0], omegaPow, gPow, gMagSqr, vPow);
        }

        if (nDimensions_ > 1)
        {
            forAllIter(List<momentFunction>, enskogFunctions_[1], iter)
            {
                (*iter)(I1s_[1], omegaPow, gPow, gMagSqr, vPow);
            }
        }
        if (nDimensions_ > 2)
        {
            forAllIter(List<momentFunction>, enskogFunctions_[2], iter)
            {
                (*iter)(I1s_[2], omegaPow, gPow, gMagSqr, vPow);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::BoltzmannCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature
)
:
    collisionKernel(dict, mesh, quadrature),
    e_(readScalar(dict.lookup("e"))),
    omega_((1.0 + e_)*0.5),
    Enskog_(dict.lookupOrDefault("Enskog", false)),
    scalarIndexes_(quadrature.nodes()[0].scalarIndexes()),
    Is_(velocityMomentOrders_.size(), velocityMomentOrders_, Zero),
    I1s_(velocityIndexes_.size()),
    Cs_(momentOrders_.size(), momentOrders_),
    gradWs_(),
    Gs_(momentOrders_.size(), momentOrders_)
{
    implicit_ = false;
    forAll(Cs_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        Cs_.set
        (
            momentOrder,
            new volScalarField
            (
                IOobject
                (
                    "collisionalSource."
                  + mappedList<scalar>::listToWord(momentOrder),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    quadrature.moments()[mi].dimensions()/dimTime,
                    Zero
                )
            )
        );
    }

    if (Enskog_)
    {
        forAll(velocityIndexes_, cmpt)
        {
            I1s_.set
            (
                cmpt,
                new mappedScalarList
                (
                    velocityMomentOrders_.size(),
                    velocityMomentOrders_,
                    Zero
                )
            );
        }

        gradWs_.resize(quadrature_.nodes().size());

        forAll(gradWs_, nodei)
        {
            gradWs_.set
            (
                nodei,
                new volVectorField
                (
                    fvc::grad(quadrature.nodes()[nodei].primaryWeight())
                )
            );
        }

        forAll(Gs_, mi)
        {
            const labelList& momentOrder = momentOrders_[mi];

            Gs_.set
            (
                mi,
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            IOobject::groupName
                            (
                                "collisionalFlux",
                                mappedScalarList::listToWord(momentOrder)
                            ),
                            quadrature_.moments()[0].group()
                        ),
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedVector
                    (
                        "zero",
                        quadrature_.moments()(momentOrder).dimensions()*dimLength/dimTime,
                        Zero
                    ),
                    wordList
                    (
                        quadrature_.moments()[0].boundaryField().size(), 
                        "zeroGradient"
                    )
                )
            );
        }
    }

    //- Check if moments exist and add coefficient functions
    mappedLabelList map(velocityMomentOrders_.size(), velocityMomentOrders_, 0);

    addIFunction1(map, 0)

    addIFunction3(map, 0,0,1)
    addIFunction2(map, 0,1)
    addIFunction1(map, 1)

    addIFunction3(map, 0,0,2)
    addIFunction3(map, 0,1,1)
    addIFunction2(map, 0,2)
    addIFunction3(map, 1,0,1)
    addIFunction2(map, 1,1)
    addIFunction1(map, 2)

    addIFunction3(map, 0,0,3)
    addIFunction3(map, 0,1,2)
    addIFunction3(map, 0,2,1)
    addIFunction2(map, 0,3)
    addIFunction3(map, 1,0,2)
    addIFunction3(map, 1,1,1)
    addIFunction2(map, 1,2)
    addIFunction3(map, 2,0,1)
    addIFunction2(map, 2,1)
    addIFunction1(map, 3)
//
    addIFunction3(map, 0,0,4)
//     addIFunction3(map, 0,1,3)
//     addIFunction3(map, 0,2,2)
//     addIFunction3(map, 0,3,1)
    addIFunction2(map, 0,4)
//     addIFunction3(map, 1,0,3)
//     addIFunction2(map, 1,3)
//     addIFunction3(map, 2,0,2)
//     addIFunction2(map, 2,2)
//     addIFunction3(map, 3,0,1)
//     addIFunction2(map, 3,1)
    addIFunction1(map, 4)
//
//     addIFunction3(map, 0,0,5)
//     addIFunction3(map, 0,1,4)
//     addIFunction3(map, 0,2,3)
//     addIFunction3(map, 0,3,2)
//     addIFunction3(map, 0,4,1)
//     addIFunction2(map, 0,5)
//     addIFunction3(map, 1,0,4)
//     addIFunction2(map, 1,4)
//     addIFunction3(map, 2,0,3)
//     addIFunction2(map, 2,3)
//     addIFunction3(map, 3,0,2)
//     addIFunction2(map, 3,2)
//     addIFunction3(map, 4,0,1)
//     addIFunction2(map, 4,1)
//     addIFunction1(map, 5)
//
//     addIFunction3(map, 0,1,5)
//     addIFunction3(map, 0,2,4)
//     addIFunction3(map, 0,4,2)
//     addIFunction3(map, 0,5,1)
//     addIFunction3(map, 1,0,5)
//     addIFunction2(map, 1,5)
//     addIFunction3(map, 2,0,4)
//     addIFunction2(map, 2,4)
//     addIFunction3(map, 4,0,2)
//     addIFunction2(map, 4,2)
//     addIFunction3(map, 5,0,1)
//     addIFunction2(map, 5,1)
//
//     addIFunction3(map, 0,2,5)
//     addIFunction3(map, 0,5,1)
//     addIFunction3(map, 2,0,5)
//     addIFunction2(map, 2,5)
//     addIFunction3(map, 5,0,2)
//     addIFunction2(map, 5,2)
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::~BoltzmannCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::
preUpdate()
{
    if (Enskog_)
    {
        forAll(gradWs_, nodei)
        {
            gradWs_[nodei] = 
                fvc::grad(quadrature_.nodes()[nodei].primaryWeight());
        }
    }
}


void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateCells
(
    const label celli
)
{
    forAll(Cs_, momenti)
    {
        Cs_[momenti][celli] = Zero;
    }

    if (Enskog_)
    {
        forAll(Gs_, momenti)
        {
            Gs_[momenti][celli] = Zero;
        }
    }

    scalar alpha = quadrature_.moments()(0)[celli];
    scalar alphac = 1.0 - alpha;

    // Monodisperse case
    if (sizeIndex_ == -1)
    {
        scalar c = min(alpha/0.63, 0.999);
        scalar g0 = (2.0 - c)/(2.0*pow3(1.0 - c)) + 1.1603*c;

        forAll(quadrature_.nodes(), nodei)
        {
            const volVelocityNode& node1 = quadrature_.nodes()[nodei];

            forAll(quadrature_.nodes(), nodej)
            {
                const volVelocityNode& node2 = quadrature_.nodes()[nodej];
                updateI(celli, nodei, nodej, omega_);

                forAll(Cs_, momenti)
                {
                    const labelList& momentOrder = momentOrders_[momenti];
                    labelList vMomentOrder(velocityIndexes_.size(), 0);

                    forAll(velocityIndexes_, cmpt)
                    {
                        vMomentOrder[cmpt] = 
                            momentOrder[velocityIndexes_[cmpt]];
                    }

                    //- Zero order source term
                    scalar cSource =
                        6.0*g0/dp_()[celli]
                       *node1.primaryWeight()[celli]
                       *node2.primaryWeight()[celli]
                       *mag
                        (
                            node1.velocityAbscissae()[celli]
                          - node2.velocityAbscissae()[celli]
                        )
                       *Is_(vMomentOrder);

                    //- Enskog term
                    if (Enskog_)
                    {
                        scalar eSource(0);

                        forAll(velocityIndexes_, m)
                        {
                            scalar I1m = I1s_[m](vMomentOrder);
                            eSource +=
                                I1m
                               *(
                                    node2.primaryWeight()[celli]
                                   *gradWs_[nodei][celli][m]
                                  - node1.primaryWeight()[celli]
                                   *gradWs_[nodej][celli][m]
                                );

                            Gs_[momenti][celli][m] +=
                                3.0*g0
                               *I1m
                               *node1.primaryWeight()[celli]
                               *node2.primaryWeight()[celli];
                        }

                        cSource += 3.0*g0*eSource;
                    }

                    // Integrate over non-velocity components
                    forAll(scalarIndexes_, cmpt)
                    {
                        scalar absCmpt =
                            pow
                            (
                                node1.primaryAbscissae()[scalarIndexes_[cmpt]][celli],
                                momentOrder[scalarIndexes_[cmpt]]
                            );

                        cSource *= absCmpt;

                        if (Enskog_)
                        {
                            Gs_[momenti][celli] *= absCmpt;
                        }
                    }

                    Cs_[momenti][celli] += cSource;
                }
            }
        }

        return;
    }

    // Polydisperse case
    forAll(quadrature_.nodes(), nodei)
    {
        const label sizei = nodeIndexes_[nodei][sizeIndex_];
        const volVelocityNode& node1 = quadrature_.nodes()[nodei];
        scalar d1 = d(sizei, celli);
        scalar V1 = Foam::constant::mathematical::pi/6.0*pow3(d1);
        scalar mass1 = V1*rhos_[sizei];
        scalar n1 = node1.primaryWeight()[celli]/V1;

        forAll(quadrature_.nodes(), nodej)
        {
            const label sizej = nodeIndexes_[nodej][sizeIndex_];
            const volVelocityNode& node2 = quadrature_.nodes()[nodej];
            scalar d2 = d(sizej, celli);
            scalar V2 = Foam::constant::mathematical::pi/6.0*pow3(d2);
            scalar mass2 = V2*rhos_[sizej];
            scalar n2 = node2.primaryWeight()[celli]/V2;

            scalar d12 = (d1 + d2)*0.5;
            scalar XiSqr = sqr(d12/d2);
            scalar omega = mass2*(1.0 + e_)/(mass1 + mass2);

            scalar xi = Foam::constant::mathematical::pi*(n1*sqr(d1) 
                + n2*sqr(d2))/6.0;

            scalar g012 =
                1.0/alphac
              + 1.5*xi*d1*d2/(sqr(alphac)*(d12))
              + 0.5*sqr(xi)/pow3(alphac)*sqr(d1*d2/d12);

            if (omega > SMALL)
            {
                updateI(celli, nodei, nodej, omega);

                forAll(Cs_, momenti)
                {
                    const labelList& momentOrder = momentOrders_[momenti];
                    labelList vMomentOrder(velocityIndexes_.size(), 0);

                    forAll(velocityIndexes_, cmpt)
                    {
                        vMomentOrder[cmpt] = 
                            momentOrder[velocityIndexes_[cmpt]];
                    }

                    //- Zero order source term
                    scalar cSource =
                        6.0*XiSqr*g012/d2
                       *node1.primaryWeight()[celli]
                       *node2.primaryWeight()[celli]
                       *mag
                        (
                            node1.velocityAbscissae()[celli]
                          - node2.velocityAbscissae()[celli]
                        )
                       *Is_(vMomentOrder);

                    //- Enskog term
                    if (Enskog_)
                    {
                        scalar enskogCoeff = 3.0*XiSqr*g012*d1/d2;
                        scalar eSource(0);
                        forAll(velocityIndexes_, m)
                        {
                            scalar I1m = I1s_[m](vMomentOrder);
                            eSource +=
                                I1m
                               *(
                                    d1*node2.primaryWeight()[celli]
                                   *gradWs_[nodei][celli][m]
                                  - d2*node1.primaryWeight()[celli]
                                   *gradWs_[nodej][celli][m]
                                );
                            Gs_[momenti][celli][m] =
                                enskogCoeff
                               *I1m
                               *node1.primaryWeight()[celli]
                               *node2.primaryWeight()[celli];
                        }
                        cSource += enskogCoeff*eSource;
                    }

                    // Integrate over non-velocity abscissae
                    forAll(scalarIndexes_, cmpt)
                    {
                        scalar absCmpt =
                            pow
                            (
                                node1.primaryAbscissae()[scalarIndexes_[cmpt]][celli],
                                momentOrder[scalarIndexes_[cmpt]]
                            );
                        cSource *= absCmpt;

                        if (Enskog_)
                        {
                            Gs_[momenti][celli] *= absCmpt;
                        }
                    }

                    Cs_[momenti][celli] += cSource;
                }
            }
        }
    }
}


Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::explicitCollisionSource
(
    const labelList& momentOrder,
    const label celli
) const
{
    return Cs_(momentOrder)[celli];
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::implicitCollisionSource(const volVelocityMoment& m) const
{
    tmp<fvScalarMatrix> iSource
    (
        new fvScalarMatrix
        (
            m,
            m.dimensions()*dimVolume/dimTime
        )
    );

    if (sizeIndex_ != -1 && Enskog_)
    {
        iSource.ref() -= fvc::div(Gs_(m.cmptOrders()))();
    }
    
    return iSource;

}


// ************************************************************************* //
