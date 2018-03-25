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
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::I1D
(
    const label celli,
    const label node1,
    const label node2
)
{
    const vector& u1 = quadrature_.nodes()[node1].primaryAbscissa()[celli];
    const vector& u2 = quadrature_.nodes()[node2].primaryAbscissa()[celli];
    scalar v1(u1.x());
    vector g = u1 - u2;
    scalar g1(g.x());
    scalar g1Sqr(sqr(g1));
    scalar gMag = mag(g);
    scalar gMagSqr = sqr(gMag);
    scalar omegaSqr = sqr(omega_);
    scalar omegaPow3 = omega_*omegaSqr;
    scalar omegaPow4 = omegaSqr*omegaSqr;

    Is_(0) = 0.0;
    Is_(1) = -(omega_/2.0)*g1;
    Is_(2) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g1Sqr - omega_*g1*v1;
    Is_(3) =
      - (omegaPow3/8.0)*(gMagSqr + g1Sqr)*g1
      + (omegaSqr/4.0)*(gMagSqr + 3.0*g1Sqr)*v1
      - (1.5*omega_)*g1*sqr(v1);
    Is_(4) =
        (omegaPow4/80.0)*(sqr(gMagSqr) + 10.0*gMagSqr*g1Sqr + 5.0*sqr(g1Sqr))
      - (omegaPow3/2.0)*(gMagSqr + g1Sqr)*g1*v1
      + (omegaSqr/2.0)*(gMagSqr + 3.0*g1Sqr)*sqr(v1)
      - 2.0*omega_*g1*pow3(v1);
}

void
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::I2D
(
    const label celli,
    const label node1,
    const label node2
)
{
    const vector& u1 = quadrature_.nodes()[node1].primaryAbscissa()[celli];
    const vector& u2 = quadrature_.nodes()[node2].primaryAbscissa()[celli];
    scalar v1(u1.x());
    scalar v2(u1.y());
    vector g = u1 - u2;
    scalar g1(g.x());
    scalar g1Sqr(sqr(g1));
    scalar g2(g.y());
    scalar g2Sqr(sqr(g2));
    scalar gMag = mag(g);
    scalar gMagSqr = sqr(gMag);
    scalar omegaSqr = sqr(omega_);
    scalar omegaPow3 = omega_*omegaSqr;
    scalar omegaPow4 = omegaSqr*omegaSqr;

    Is_(0,0) = 0.0;
    Is_(1,0) = -(omega_/2.0)*g1;
    Is_(0,1) = -(omega_/2.0)*g2;
    Is_(2,0) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g1Sqr - omega_*g1*v1;
    Is_(1,1) = (omegaSqr/4.0)*g1*g2 - (omega_/2.0)*(v1*g2 + g1*v2);
    Is_(0,2) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g2Sqr - omega_*g2*v2;
    Is_(3,0) =
      - (omegaPow3/8.0)*(gMagSqr + g1Sqr)*g1
      + (omegaSqr/4.0)*(gMagSqr + 3.0*g1Sqr)*v1
      - (1.5*omega_)*g1*sqr(v1);
    Is_(0,3) =
      - (omegaPow3/8.0)*(gMagSqr + g2Sqr)*g2
      + (omegaSqr/4.0)*(gMagSqr + 3.0*g2Sqr)*v2
      - (1.5*omega_)*g2*sqr(v2);
    Is_(4,0) =
        (omegaPow4/80.0)*(sqr(gMagSqr) + 10.0*gMagSqr*g1Sqr + 5.0*sqr(g1Sqr))
      - (omegaPow3/2.0)*(gMagSqr + g1Sqr)*g1*v1
      + (omegaSqr/2.0)*(gMagSqr + 3.0*g1Sqr)*sqr(v1)
      - 2.0*omega_*g1*pow3(v1);
    Is_(0,4) =
        (omegaPow4/80.0)*(sqr(gMagSqr) + 10.0*gMagSqr*g2Sqr + 5.0*sqr(g2Sqr))
      - (omegaPow3/2.0)*(gMagSqr + g2Sqr)*g2*v2
      + (omegaSqr/2.0)*(gMagSqr + 3.0*g2Sqr)*sqr(v2)
      - 2.0*omega_*g2*pow3(v2);
}

void
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::I3D
(
    const label celli,
    const label node1,
    const label node2
)
{
    const vector& u1 = quadrature_.nodes()[node1].primaryAbscissa()[celli];
    const vector& u2 = quadrature_.nodes()[node2].primaryAbscissa()[celli];
    scalar v1(u1.x());
    scalar v2(u1.y());
    scalar v3(u1.z());
    vector g = u1 - u2;
    scalar g1(g.x());
    scalar g1Sqr(sqr(g1));
    scalar g2(g.y());
    scalar g2Sqr(sqr(g2));
    scalar g3(g.z());
    scalar g3Sqr(sqr(g3));
    scalar gMag = mag(g);
    scalar gMagSqr = sqr(gMag);
    scalar omegaSqr = sqr(omega_);
    scalar omegaPow3 = omega_*omegaSqr;
    scalar omegaPow4 = omegaSqr*omegaSqr;

    Is_(0,0,0) = 0.0;
    Is_(1,0,0) = -(omega_/2.0)*g1;
    Is_(0,1,0) = -(omega_/2.0)*g2;
    Is_(0,0,1) = -(omega_/2.0)*g3;
    Is_(2,0,0) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g1Sqr - omega_*g1*v1;
    Is_(1,1,0) = (omegaSqr/4.0)*g1*g2 - (omega_/2.0)*(v1*g2 + g1*v2);
    Is_(1,0,1) = (omegaSqr/4.0)*g1*g3 - (omega_/2.0)*(v1*g3 + g1*v3);
    Is_(0,2,0) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g2Sqr - omega_*g2*v2;
    Is_(0,1,1) = (omegaSqr/4.0)*g2*g3 - (omega_/2.0)*(v2*g3 + g2*v3);
    Is_(0,0,2) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g3Sqr - omega_*g3*v3;
    Is_(3,0,0) =
      - (omegaPow3/8.0)*(gMagSqr + g1Sqr)*g1
      + (omegaSqr/4.0)*(gMagSqr + 3.0*g1Sqr)*v1
      - (1.5*omega_)*g1*sqr(v1);
    Is_(0,3,0) =
      - (omegaPow3/8.0)*(gMagSqr + g2Sqr)*g2
      + (omegaSqr/4.0)*(gMagSqr + 3.0*g2Sqr)*v2
      - (1.5*omega_)*g2*sqr(v2);
    Is_(0,0,3) =
      - (omegaPow3/8.0)*(gMagSqr + g3Sqr)*g3
      + (omegaSqr/4.0)*(gMagSqr + 3.0*g3Sqr)*v3
      - (1.5*omega_)*g3*sqr(v3);
    Is_(4,0,0) =
        (omegaPow4/80.0)*(sqr(gMagSqr) + 10.0*gMagSqr*g1Sqr + 5.0*sqr(g1Sqr))
      - (omegaPow3/2.0)*(gMagSqr + g1Sqr)*g1*v1
      + (omegaSqr/2.0)*(gMagSqr + 3.0*g1Sqr)*sqr(v1)
      - 2.0*omega_*g1*pow3(v1);
    Is_(0,4,0) =
        (omegaPow4/80.0)*(sqr(gMagSqr) + 10.0*gMagSqr*g2Sqr + 5.0*sqr(g2Sqr))
      - (omegaPow3/2.0)*(gMagSqr + g2Sqr)*g2*v2
      + (omegaSqr/2.0)*(gMagSqr + 3.0*g2Sqr)*sqr(v2)
      - 2.0*omega_*g2*pow3(v2);
    Is_(0,0,4) =
        (omegaPow4/80.0)*(sqr(gMagSqr) + 10.0*gMagSqr*g3Sqr + 5.0*sqr(g3Sqr))
      - (omegaPow3/2.0)*(gMagSqr + g3Sqr)*g3*v3
      + (omegaSqr/2.0)*(gMagSqr + 3.0*g3Sqr)*sqr(v3)
      - 2.0*omega_*g3*pow3(v3);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::BoltzmannCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature,
    const bool ode
)
:
    collisionKernel(dict, mesh, quadrature, ode),
    dp_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                "d",
                quadrature.moments()[0].group()
            )
        )
    ),
    e_(dict.lookupType<scalar>("e")),
    omega_((1.0 + e_)*0.5),
    Is_(momentOrders_.size(), momentOrders_),
    Cs_(momentOrders_.size(), momentOrders_)
{
    Info<<"found d.particles"<<endl;
    if (!ode)
    {
        FatalErrorInFunction
            << "Boltzmann collision kernel does not support implicit" << nl
            << "solutions to the collisional source term."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::~BoltzmannCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateCells
(
    const label celli
)
{
    forAll(Cs_, momenti)
    {
        Cs_[momenti] = 0.0;
    }

    forAll(quadrature_.nodes(), nodei)
    {
        const volVectorNode& node1 = quadrature_.nodes()[nodei];
        forAll(quadrature_.nodes(), nodej)
        {
            const volVectorNode& node2 = quadrature_.nodes()[nodej];

            if (nDimensions_ == 1)
            {
                I1D(celli, nodei, nodej);
            }
            else if (nDimensions_ == 2)
            {
                I2D(celli, nodei, nodej);
            }
            else if (nDimensions_ == 3)
            {
                I3D(celli, nodei, nodej);
            }

            forAll(Cs_, momenti)
            {
                const labelList& momentOrder = momentOrders_[momenti];
                Cs_(momentOrder) +=
                    node1.primaryWeight()[celli]
                   *node2.primaryWeight()[celli]
                   *mag
                    (
                        node1.primaryAbscissa()[celli]
                      - node2.primaryAbscissa()[celli]
                    )
                   *Is_(momentOrder);
            }
        }
    }

    scalar alpha = quadrature_.moments()(labelList(nDimensions_, 0))[celli];
    scalar g0 = 1.0/(1 - alpha)
      + 3*alpha/(2*sqr(1 - alpha))
      + sqr(alpha)/(2*pow3(1 - alpha));

    forAll(Cs_, momenti)
    {
        Cs_[momenti] *= 6.0*g0/dp_[celli];
    }
}

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateFields()
{
    NotImplemented;
}

Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::explicitCollisionSource(const label mi, const label celli) const
{
    return Cs_[mi];
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::implicitCollisionSource(const volVectorMoment& m) const
{
    NotImplemented;

    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            m,
            m.dimensions()*dimVolume/dimTime
        )
    );
}
// ************************************************************************* //
