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

\*---------------------------------------------------------------------------*/

#include "BoltzmannCollision.H"
#include "hyperbolicConditionalMomentInversion.H"
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

    // Store powers of velocities
    vector g = u1 - u2;
    scalarList omegaPow(6, omega);
    vectorList gPow(6, g);
    scalarList gMagPow(6, mag(g));
    vectorList vPow(6, u1);

    for (label powi = 2; powi < gPow.size(); powi++)
    {
        omegaPow[powi] = pow(omegaPow[powi], powi);
        gMagPow[powi] = pow(gMagPow[powi], powi);
        for (label cmpt = 0; cmpt < nDimensions_; cmpt++)
        {
            gPow[powi][cmpt] = pow(gPow[powi][cmpt], powi);
            vPow[powi][cmpt] = pow(vPow[powi][cmpt], powi);
        }
    }

    forAllIter(List<momentFunction>, coefficientFunctions_, iter)
    {
        (*iter)(Is_, omegaPow, gPow, gMagPow, vPow);
    }

    if (Enskog_)
    {
        forAllIter(List<momentFunction>, enskogFunctions_[0], iter)
        {
            (*iter)(I1s_[0], omegaPow, gPow, gMagPow, vPow);
        }

        if (nDimensions_ > 1)
        {
            forAllIter(List<momentFunction>, enskogFunctions_[1], iter)
            {
                (*iter)(I1s_[1], omegaPow, gPow, gMagPow, vPow);
            }
        }
        if (nDimensions_ > 2)
        {
            forAllIter(List<momentFunction>, enskogFunctions_[2], iter)
            {
                (*iter)(I1s_[2], omegaPow, gPow, gMagPow, vPow);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::BoltzmannCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature
)
:
    collisionKernel(dict, mesh, quadrature),
    e_(dict.lookupType<scalar>("e")),
    omega_((1.0 + e_)*0.5),
    Enskog_(dict.lookupOrDefault("Enskog", false)),
    scalarIndexes_(quadrature.nodes()[0].scalarIndexes()),
    Is_(velocityMomentOrders_.size(), velocityMomentOrders_, 0.0),
    I1s_(velocityIndexes_.size()),
    Cs_(momentOrders_.size(), momentOrders_),
    dp_
    (
        nSizes_ <= 0
      ? lookupOrInitialize
        (
            mesh,
            IOobject::groupName("d", quadrature.moments()[0].group()),
            dict,
            "d",
            dimLength
        )
      : tmp<volScalarField>()
    ),
    rhop_
    (
        lookupOrInitialize
        (
            mesh,
            IOobject::groupName("thermo:rho", quadrature.moments()[0].group()),
            dict,
            "rho",
            dimDensity
        )
    ),
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
                    0.0
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
                    0.0
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
                    wordList(quadrature_.moments()[0].boundaryField().size(), "zeroGradient")
                )
            );
        }
    }

    mappedLabelList map(velocityMomentOrders_.size(), velocityMomentOrders_, 0);
//
    addIFunction1(map, 0)
//
    addIFunction3(map, 0,0,1)
    addIFunction2(map, 0,1)
    addIFunction1(map, 1)
//
    addIFunction3(map, 0,0,2)
    addIFunction3(map, 0,1,1)
    addIFunction2(map, 0,2)
    addIFunction3(map, 1,0,1)
    addIFunction2(map, 1,1)
    addIFunction1(map, 2)
//
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
            gradWs_[nodei] = fvc::grad(quadrature_.nodes()[nodei].primaryWeight());
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
        Cs_[momenti][celli] = 0.0;
    }
    if (Enskog_)
    {
        forAll(Gs_, momenti)
        {
            Gs_[momenti][celli] = Zero;
        }
    }

    scalar alpha = quadrature_.moments()(labelList(nDimensions_, 0))[celli];
    scalar c = min(alpha/0.63, 0.999);
    scalar alphac = 1.0 - alpha;
    scalar g0 = (2.0 - c)/(2.0*pow3(1.0 - c));

    if (sizeIndex_ == -1)
    {
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
                        vMomentOrder[cmpt] = momentOrder[velocityIndexes_[cmpt]];
                    }

                    //- Zero order
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
                        scalar eSource = 0.0;
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

    scalar alphard = 0.0;

    forAll(quadrature_.nodes(), nodei)
    {
        const volVelocityNode& node = quadrature_.nodes()[nodei];
        alphard +=
            node.primaryWeight()[celli]
            /max(node.primaryAbscissae()[sizeIndex_][celli], 1e-6);
    }

    scalar pi = Foam::constant::mathematical::pi;

    forAll(quadrature_.nodes(), nodei)
    {
        const volVelocityNode& node1 = quadrature_.nodes()[nodei];
        forAll(quadrature_.nodes(), nodej)
        {
            const volVelocityNode& node2 = quadrature_.nodes()[nodej];


            scalar d1 = max(node1.primaryAbscissae()[sizeIndex_][celli], 1e-10);
            scalar d2 = max(node2.primaryAbscissae()[sizeIndex_][celli], 1e-10);
            scalar d12 = (d1 + d2)*0.5;
            scalar XiSqr = sqr(d12/d2);
            scalar mass1 = pi/6.0*pow3(d1)*rhop_()[celli];
            scalar mass2 = pi/6.0*pow3(d2)*rhop_()[celli];
            scalar omega = mass2*(1.0 + e_)/max(mass1 + mass2, small);
            scalar g012 = 1.0/alphac + 3.0*d1*d2*alphard/(sqr(alphac)*(d1 + d2));

            if (omega > 1e-10)
            {
                updateI(celli, nodei, nodej, omega);

                forAll(Cs_, momenti)
                {
                    const labelList& momentOrder = momentOrders_[momenti];
                    labelList vMomentOrder(velocityIndexes_.size(), 0);
                    forAll(velocityIndexes_, cmpt)
                    {
                        vMomentOrder[cmpt] = momentOrder[velocityIndexes_[cmpt]];
                    }

                    //- Zero order
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
                        scalar eSource = 0.0;
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
