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

    scalar v1(u1.x());
    scalar v2(u1.y());
    scalar v3(u1.z());
    vector g = u1 - u2;
    scalar g1(g.x());
    scalar g1Sqr(sqr(g1));
    scalar g1Pow4(sqr(g1Sqr));
    scalar g2(g.y());
    scalar g2Sqr(sqr(g2));
    scalar g2Pow4(sqr(g2Sqr));
    scalar g3(g.z());
    scalar g3Sqr(sqr(g3));
    scalar g3Pow4(sqr(g3Sqr));
    scalar gMag = mag(g);
    scalar gMagSqr = sqr(gMag);
    scalar gPow4 = sqr(gMagSqr);
    scalar omegaSqr = sqr(omega);
    scalar omegaPow3 = omega*omegaSqr;
    scalar omegaPow4 = omegaSqr*omegaSqr;

    Is_(1) = -(omega/2.0)*g1;
    Is_(2) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g1Sqr - omega*g1*v1;
    Is_(3) =
      - (omegaPow3/8.0)*(gMagSqr + g1Sqr)*g1
      + (omegaSqr/4.0)*(gMagSqr + 3.0*g1Sqr)*v1
      - (1.5*omega)*g1*sqr(v1);
    Is_(4) =
        (omegaPow4/80.0)*(gPow4 + 10.0*gMagSqr*g1Sqr + 5.0*g1Pow4)
      - (omegaPow3/2.0)*(gMagSqr + g1Sqr)*g1*v1
      + (omegaSqr/2.0)*(gMagSqr + 3.0*g1Sqr)*sqr(v1)
      - 2.0*omega*g1*pow3(v1);

    if (nDimensions_ > 1)
    {
        Is_(0,1) = -(omega/2.0)*g2;
        Is_(1,1) = (omegaSqr/4.0)*g1*g2 - (omega/2.0)*(v1*g2 + g1*v2);
        Is_(0,2) =
            (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g2Sqr - omega*g2*v2;
        Is_(0,3) =
          - (omegaPow3/8.0)*(gMagSqr + g2Sqr)*g2
          + (omegaSqr/4.0)*(gMagSqr + 3.0*g2Sqr)*v2
          - (1.5*omega)*g2*sqr(v2);
        Is_(0,4) =
            (omegaPow4/80.0)*(gPow4 + 10.0*gMagSqr*g2Sqr + 5.0*g2Pow4)
          - (omegaPow3/2.0)*(gMagSqr + g2Sqr)*g2*v2
          + (omegaSqr/2.0)*(gMagSqr + 3.0*g2Sqr)*sqr(v2)
          - 2.0*omega*g2*pow3(v2);
    }
    if (nDimensions_ > 2)
    {
        Is_(0,0,1) = -(omega/2.0)*g3;
        Is_(1,0,1) = (omegaSqr/4.0)*g1*g3 - (omega/2.0)*(v1*g3 + g1*v3);
        Is_(0,1,1) = (omegaSqr/4.0)*g2*g3 - (omega/2.0)*(v2*g3 + g2*v3);
        Is_(0,0,2) =
            (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g3Sqr - omega*g3*v3;
        Is_(0,0,3) =
          - (omegaPow3/8.0)*(gMagSqr + g3Sqr)*g3
          + (omegaSqr/4.0)*(gMagSqr + 3.0*g3Sqr)*v3
          - (1.5*omega)*g3*sqr(v3);
        Is_(0,0,4) =
            (omegaPow4/80.0)*(gPow4 + 10.0*gMagSqr*g3Sqr + 5.0*g3Pow4)
          - (omegaPow3/2.0)*(gMagSqr + g3Sqr)*g3*v3
          + (omegaSqr/2.0)*(gMagSqr + 3.0*g3Sqr)*sqr(v3)
          - 2.0*omega*g3*pow3(v3);
    }

    if (Enskog_)
    {
        I1s_[0](1) = (2.0*omega/15.0)*(gMagSqr + 2.0*g1Sqr);
        I1s_[0](2) =
          - (2.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g1Sqr)*g1
          + (4.0*omega/15.0)*(gMagSqr + 2.0*g1Sqr)*v1;
        I1s_[0](3) =
            (2.0*omegaSqr/315.0)
           *(3.0*gPow4 + 24.0*gMagSqr*g1Sqr + 8.0*pow4(g1))
          - (6.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g1Sqr)*g1*v1
          + (2.0*omega/5.0)*(gMagSqr + 2.0*g1Sqr)*sqr(v1);
        I1s_[0](4) =
          - (2.0*omegaPow4/693.0)
           *(15.0*gPow4 + 40.0*gMagSqr*g1Sqr + 8.0*pow4(g1))*g1
          + (8.0*omegaPow3/315.0)*(3.0*gPow4 + 24.0*gMagSqr*g1Sqr + 8.0*g1Pow4)*v1
          - (12.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g1Sqr)*g1*sqr(v1)
          + (8.0*omega/15.0)*(gMagSqr + 2.0*g1Sqr)*pow3(v1);

        if (nDimensions_ > 1)
        {
            I1s_[1](1) = (4.0*omega/15.0)*g1*g2;

            I1s_[0](0, 1) = (4.0*omega/15.0)*g1*g2;
            I1s_[1](0, 1) = (2.0*omega/15.0)*(gMagSqr + 2.0*g2Sqr);

            I1s_[0](1,1) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g2
              + (4.0*omega/15.0)*g1*g2*v1
              + (2.0*omega/15.0)*(gMagSqr + 2.0*g1Sqr)*v2;
            I1s_[1](1,1) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g1
              + (4.0*omega/15.0)*g2*g1*v2
              + (2.0*omega/15.0)*(gMagSqr + 2.0*g2Sqr)*v1;

            I1s_[1](2) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g2
              + (8.0*omega/15.0)*g2*g1*v1;

            I1s_[0](0,2) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g1
              + (8.0*omega/15.0)*g1*g2*v2;
            I1s_[1](0,2) =
              - (2.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g2Sqr)*g2
              + (4.0*omega/15.0)*(gMagSqr + 2.0*g2Sqr)*v2;

            I1s_[1](3) =
                (8.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g1Sqr)*g2*g1
              - (6.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g2*v1
              + (4.0*omega/5.0)*g1*g2*sqr(v1);

            I1s_[0](0,3) =
                (8.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g2Sqr)*g1*g2
              - (6.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g1*v2
              + (4.0*omega/5.0)*g1*g2*sqr(v2);
            I1s_[1](0,3) =
                (2.0*omegaSqr/315.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g2Sqr + 8.0*g2Pow4)
              - (6.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g2Sqr)*g2*v2
              + (2.0*omega/5.0)*(gMagSqr + 2.0*g2Sqr)*sqr(v2);

            I1s_[1](4) =
              - (2.0*omegaPow4/693.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g1Sqr + 8.0*g1Pow4)*g2
              + (32.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g1Sqr)*g2*g1*v1
              - (12.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g2*sqr(v1)
              + (16.0*omega/15.0)*g2*g1*pow3(v1);

            I1s_[0](0,4) =
              - (2.0*omegaPow4/693.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g2Sqr + 8.0*g2Pow4)*g1
              + (32.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g2Sqr)*g1*g2*v2
              - (12.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g1*sqr(v2)
              + (16.0*omega/15.0)*g1*g2*pow3(v2);
            I1s_[1](0,4) =
              - (2.0*omegaPow4/693.0)
               *(15.0*gPow4 + 40.0*gMagSqr*g2Sqr + 8.0*g2Pow4)*g2
              + (8.0*omegaPow3/315.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g2Sqr + 8.0*g2Pow4)*v2
              - (12.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g2Sqr)*g2*sqr(v2)
              + (8.0*omega/15.0)*(gMagSqr + 2.0*g2Sqr)*pow3(v2);
        }
        if (nDimensions_ > 2)
        {
            I1s_[2](1) = (4.0*omega/15.0)*g1*g3;
            I1s_[2](0,1) = (4.0*omega/15.0)*g2*g3;

            I1s_[0](0,0,1) = (4.0*omega/15.0)*g1*g3;
            I1s_[1](0,0,1) = (4.0*omega/15.0)*g2*g3;
            I1s_[2](0,0,1) = (2.0*omega/15.0)*(gMagSqr + 2.0*g3Sqr);

            I1s_[2](1,1) =
              - (4.0*omegaSqr/35.0)*g1*g2*g3
              + (4.0*omega/15.0)*g3*(g1*v2 + g2*v1);


            I1s_[0](1,0,1) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g3
              + (4.0*omega/15.0)*g1*g3*v1
              + (2.0*omega/15.0)*(gMagSqr + 2.0*g1Sqr)*v3;
            I1s_[1](1,0,1) =
              - (4.0*omegaSqr/35.0)*g1*g2*g3
              + (4.0*omega/15.0)*g2*(g1*v3 + g3*v1);
            I1s_[2](1,0,1) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g1
              + (4.0*omega/15.0)*g3*g1*v3
              + (2.0*omega/15.0)*(gMagSqr + 2.0*g3Sqr)*v1;

            I1s_[0](0,1,1) =
              - (4.0*omegaSqr/35.0)*g1*g2*g3
              + (4.0*omega/15.0)*g1*(g2*v3 + g3*v2);
            I1s_[1](0,1,1) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g3
              + (4.0*omega/15.0)*g2*g3*v2
              + (2.0*omega/15.0)*(gMagSqr + 2.0*g2Sqr)*v3;
            I1s_[2](0,1,1) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g2
              + (4.0*omega/15.0)*g2*g3*v3
              + (2.0*omega/15.0)*(gMagSqr + 2.0*g3Sqr)*v2;

            I1s_[2](2) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g3
              + (8.0*omega/15.0)*g1*g3*v1;

            I1s_[2](0, 2) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g3
              + (8.0*omega/15.0)*g2*g3*v2;

            I1s_[0](0,0,2) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g1
              + (8.0*omega/15.0)*g1*g3*v3;
            I1s_[1](0,0,2) =
              - (2.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g2
              + (8.0*omega/15.0)*g2*g3*v3;
            I1s_[2](0,0,2) =
              - (2.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g3Sqr)*g3
              + (4.0*omega/15.0)*(gMagSqr + 2.0*g3Sqr)*v3;

            I1s_[2](3) =
                (8.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g1Sqr)*g3*g1
              - (6.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g3*v1
              + (4.0*omega/5.0)*g1*g3*sqr(v1);

            I1s_[2](0,3) =
                (8.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g2Sqr)*g3*g2
              - (6.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g3*v2
              + (4.0*omega/5.0)*g3*g2*sqr(v2);

            I1s_[0](0,0,3) =
                (8.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g3Sqr)*g1*g3
              - (6.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g1*v3
              + (4.0*omega/5.0)*g1*g3*sqr(v3);
            I1s_[1](0,0,3) =
                (8.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g3Sqr)*g2*g3
              - (6.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g2*v3
              + (4.0*omega/5.0)*g2*g3*sqr(v3);
            I1s_[2](0,0,3) =
                (2.0*omegaSqr/315.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g3Sqr + 8.0*pow4(g3))
              - (6.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g3Sqr)*g3*v3
              + (2.0*omega/5.0)*(gMagSqr + 2.0*g3Sqr)*sqr(v3);

            I1s_[2](4) =
              - (2.0*omegaPow4/693.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g1Sqr + 8.0*g1Pow4)*g3
              + (32.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g1Sqr)*g3*g1*v1
              - (12.0*omegaSqr/35.0)*(gMagSqr + 2.0*g1Sqr)*g3*sqr(v1)
              + (16.0*omega/15.0)*g3*g1*pow3(v1);

            I1s_[2](0,4) =
              - (2.0*omegaPow4/693.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g2Sqr + 8.0*g2Pow4)*g3
              + (32.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g2Sqr)*g3*g2*v2
              - (12.0*omegaSqr/35.0)*(gMagSqr + 2.0*g2Sqr)*g3*sqr(v2)
              + (16.0*omega/15.0)*g3*g2*pow3(v2);


            I1s_[0](0,0,4) =
              - (2.0*omegaPow4/693.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g3Sqr + 8.0*g3Pow4)*g1
              + (32.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g3Sqr)*g1*g3*v3
              - (12.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g1*sqr(v3)
              + (16.0*omega/15.0)*g1*g3*pow3(v3);
            I1s_[1](0,0,4) =
              - (2.0*omegaPow4/693.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g3Sqr + 8.0*g3Pow4)*g2
              + (32.0*omegaPow3/315.0)*(3.0*gMagSqr + 2.0*g3Sqr)*g2*g3*v3
              - (12.0*omegaSqr/35.0)*(gMagSqr + 2.0*g3Sqr)*g2*sqr(v3)
              + (16.0*omega/15.0)*g2*g3*pow3(v3);
            I1s_[2](0,0,4) =
              - (2.0*omegaPow4/693.0)
               *(15.0*gPow4 + 40.0*gMagSqr*g3Sqr + 8.0*g3Pow4)*g3
              + (8.0*omegaPow3/315.0)
               *(3.0*gPow4 + 24.0*gMagSqr*g3Sqr + 8.0*g3Pow4)*v3
              - (12.0*omegaSqr/35.0)*(3.0*gMagSqr + 2.0*g3Sqr)*g3*sqr(v3)
              + (8.0*omega/15.0)*(gMagSqr + 2.0*g3Sqr)*pow3(v3);
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
    sizeIndex_(quadrature.nodes()[0].sizeIndex()),
    velocityIndexes_(quadrature.nodes()[0].velocityIndexes()),
    velocityMomentOrders_
    (
        multivariateMomentInversions::CHyQMOM::getMomentOrders
        (
            velocityIndexes_.size()
        )
    ),
    Is_(velocityMomentOrders_.size(), velocityMomentOrders_, 0.0),
    I1s_(velocityIndexes_.size()),
    Cs_(momentOrders_.size(), momentOrders_),
    dp_
    (
        lookupOrInitialize
        (
            mesh,
            IOobject::groupName("d", quadrature.moments()[0].group()),
            dict,
            "d",
            dimLength
        )
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
