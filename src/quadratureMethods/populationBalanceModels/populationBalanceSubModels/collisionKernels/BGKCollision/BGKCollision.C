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

#include "BGKCollision.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{
    defineTypeNameAndDebug(BGKCollision, 0);

    addToRunTimeSelectionTable
    (
        collisionKernel,
        BGKCollision,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateCells1D(const label celli)
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();

    scalar m0 = max(moments(0)[celli], SMALL);

    // Mean velocity
    scalar u = meanVelocity(m0, moments(1)[celli]);
    scalar uSqr = sqr(u);

    // Variances of velocities
    scalar sigma = max(moments(2)[celli]/m0 - uSqr, 0.0);

    Meq_(0) = moments(0)[celli];
    Meq_(1) = moments(1)[celli];
    Meq_(2) = moments(2)[celli];
    Meq_(3) = m0*(3.0*sigma*u + u*uSqr);
    Meq_(4) = m0*(6.0*uSqr*sigma + 3.0*sqr(sigma) + uSqr*uSqr);
}

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateCells2D(const label celli)
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    scalar m00 = max(moments(0,0)[celli], SMALL);

    // Mean velocity
    scalar u = meanVelocity(m00, moments(1,0)[celli]);
    scalar v = meanVelocity(m00, moments(0,1)[celli]);
    scalar uSqr = sqr(u);
    scalar vSqr = sqr(v);

    // Variances of velocities
    scalar sigma11 = max(moments(2,0)[celli]/m00 - uSqr, 0.0);
    scalar sigma22 = max(moments(0,2)[celli]/m00 - vSqr, 0.0);

    Meq_(0,0) = moments(0,0)[celli];
    Meq_(1,0) = moments(1,0)[celli];
    Meq_(0,1) = moments(0,1)[celli];
    Meq_(2,0) = moments(0,2)[celli];
    Meq_(1,1) = moments(1,1)[celli];
    Meq_(0,2) = moments(0,2)[celli];
    Meq_(3,0) = m00*(3.0*sigma11*u + u*uSqr);
    Meq_(0,3) = m00*(3.0*sigma22*v + v*vSqr);
    Meq_(4,0) = m00*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);
    Meq_(0,4) = m00*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
}

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateCells3D(const label celli)
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    scalar m000 = max(moments(0,0,0)[celli], SMALL);

    // Mean velocity
    scalar u = meanVelocity(m000, moments(1,0,0)[celli]);
    scalar v = meanVelocity(m000, moments(0,1,0)[celli]);
    scalar w = meanVelocity(m000, moments(0,0,1)[celli]);
    scalar uSqr = sqr(u);
    scalar vSqr = sqr(v);
    scalar wSqr = sqr(w);

    // Variances of velocities
    scalar sigma11 = max(moments(2,0,0)[celli]/m000 - uSqr, 0.0);
    scalar sigma22 = max(moments(0,2,0)[celli]/m000 - vSqr, 0.0);
    scalar sigma33 = max(moments(0,0,2)[celli]/m000 - wSqr, 0.0);

    Meq_(0,0,0) = moments(0,0,0)[celli];
    Meq_(1,0,0) = moments(1,0,0)[celli];
    Meq_(0,1,0) = moments(0,1,0)[celli];
    Meq_(0,0,1) = moments(0,0,1)[celli];
    Meq_(2,0,0) = moments(2,0,0)[celli];
    Meq_(1,1,0) = moments(1,1,0)[celli];
    Meq_(1,0,1) = moments(1,0,1)[celli];
    Meq_(0,2,0) = moments(0,2,0)[celli];
    Meq_(0,1,1) = moments(0,1,1)[celli];
    Meq_(0,0,2) = moments(0,0,2)[celli];
    Meq_(3,0,0) = m000*(3.0*sigma11*u + u*uSqr);
    Meq_(0,3,0) = m000*(3.0*sigma22*v + v*vSqr);
    Meq_(0,0,3) = m000*(3.0*sigma33*w + w*wSqr);
    Meq_(4,0,0) = m000*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);
    Meq_(0,4,0) = m000*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
    Meq_(0,0,4) = m000*(6.0*wSqr*sigma33 + 3.0*sqr(sigma33) + wSqr*wSqr);
}


void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateFields1D()
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    volScalarField m0(max(moments(0), SMALL));

    // Mean velocity
    volScalarField u(meanVelocity(m0, moments(1)));
    volScalarField uSqr(sqr(u));

    // Variances of velocities
    dimensionedScalar zeroVar("zero", sqr(dimVelocity), 0.0);
    volScalarField sigma(max(moments(2)/m0 - uSqr, zeroVar));

    Meqf_(0) = moments(0);
    Meqf_(1) = moments(1);
    Meqf_(2) = moments(2);
    Meqf_(3) = m0*(3.0*sigma*u + u*uSqr);
    Meqf_(4) = m0*(6.0*uSqr*sigma + 3.0*sqr(sigma) + uSqr*uSqr);
}

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateFields2D()
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    volScalarField m00(max(moments(0,0), SMALL));

    // Mean velocity
    volScalarField u(meanVelocity(m00, moments(1,0)));
    volScalarField v(meanVelocity(m00, moments(0,1)));
    volScalarField uSqr(sqr(u));
    volScalarField vSqr(sqr(v));

    // Variances of velocities
    dimensionedScalar zeroVar("zero", sqr(dimVelocity), 0.0);
    volScalarField sigma11(max(moments(2,0)/m00 - uSqr, zeroVar));
    volScalarField sigma22(max(moments(0,2)/m00 - vSqr, zeroVar));

    Meqf_(0,0) = moments(0,0);
    Meqf_(1,0) = moments(1,0);
    Meqf_(0,1) = moments(0,1);
    Meqf_(2,0) = moments(2,0);
    Meqf_(1,1) = moments(1,1);
    Meqf_(0,2) = moments(0,2);
    Meqf_(3,0) = m00*(3.0*sigma11*u + u*uSqr);
    Meqf_(0,3) = m00*(3.0*sigma22*v + v*vSqr);
    Meqf_(4,0) = m00*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);
    Meqf_(0,4) = m00*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
}

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateFields3D()
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    volScalarField m000(max(moments(0,0,0), SMALL));

    // Mean velocity
    volScalarField u(meanVelocity(m000, moments(1,0,0)));
    volScalarField v(meanVelocity(m000, moments(0,1,0)));
    volScalarField w(meanVelocity(m000, moments(0,0,1)));
    volScalarField uSqr(sqr(u));
    volScalarField vSqr(sqr(v));
    volScalarField wSqr(sqr(w));

    // Variances of velocities
    dimensionedScalar zeroVar("zero", sqr(dimVelocity), 0.0);
    volScalarField sigma11(max(moments(2,0,0)/m000 - uSqr, zeroVar));
    volScalarField sigma22(max(moments(0,2,0)/m000 - vSqr, zeroVar));
    volScalarField sigma33(max(moments(0,0,2)/m000 - wSqr, zeroVar));

    Meqf_(0,0,0) = moments(0,0,0);
    Meqf_(1,0,0) = moments(1,0,0);
    Meqf_(0,1,0) = moments(0,1,0);
    Meqf_(0,0,1) = moments(0,0,1);
    Meqf_(2,0,0) = moments(2,0,0);
    Meqf_(1,1,0) = moments(1,1,0);
    Meqf_(1,0,1) = moments(1,0,1);
    Meqf_(0,2,0) = moments(0,2,0);
    Meqf_(0,1,1) = moments(0,1,1);
    Meqf_(0,0,2) = moments(0,0,2);
    Meqf_(3,0,0) = m000*(3.0*sigma11*u + u*uSqr);
    Meqf_(0,3,0) = m000*(3.0*sigma22*v + v*vSqr);
    Meqf_(0,0,3) = m000*(3.0*sigma33*w + w*wSqr);
    Meqf_(4,0,0) = m000*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);
    Meqf_(0,4,0) = m000*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
    Meqf_(0,0,4) = m000*(6.0*wSqr*sigma33 + 3.0*sqr(sigma33) + wSqr*wSqr);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BGKCollision::BGKCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature,
    const bool ode
)
:
    collisionKernel(dict, mesh, quadrature, ode),
    tauCollisional_(dict.lookup("tau")),
    Meqf_(quadrature.moments().size(), momentOrders_),
    Meq_(quadrature.moments().size(), momentOrders_)
{
    if (!ode)
    {
//         Meqf_.setSize(quadrature.moments().size());

        forAll(Meqf_, mi)
        {
            const labelList& momentOrder = momentOrders_[mi];
            Meqf_.set
            (
                momentOrder,
                new volScalarField
                (
                    IOobject
                    (
                        "Meq" + mappedList<scalar>::listToWord(momentOrder),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "zero",
                        quadrature.moments()[mi].dimensions(),
                        0.0
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::~BGKCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateCells
(
    const label celli
)
{
    if (nDimensions_ == 1)
    {
        updateCells1D(celli);
    }
    else if (nDimensions_ == 2)
    {
        updateCells2D(celli);
    }
    else if (nDimensions_ == 3)
    {
        updateCells3D(celli);
    }
}

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateFields()
{
    if (nDimensions_ == 1)
    {
        updateFields1D();
    }
    else if (nDimensions_ == 2)
    {
        updateFields2D();
    }
    else if (nDimensions_ == 3)
    {
        updateFields3D();
    }
}

Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::explicitCollisionSource(const label mi, const label celli) const
{
    return (quadrature_.moments()[mi][celli] - Meq_[mi])/tauCollisional_.value();
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::implicitCollisionSource(const volVectorMoment& m) const
{
    return
    (
        Meqf_(m.cmptOrders())/tauCollisional_
      - fvm::Sp(1/tauCollisional_, m)
    );
}
// ************************************************************************* //
