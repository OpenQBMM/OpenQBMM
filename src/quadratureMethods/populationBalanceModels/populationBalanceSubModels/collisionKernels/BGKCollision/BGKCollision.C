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
::updateCells(const label celli)
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    scalar m0 = max(moments(0)[celli], SMALL);

    // Mean velocity
    scalar u = meanVelocity(m0, moments(1)[celli]);
    scalar uSqr = sqr(u);
    scalar sigma11 = max(moments(2)[celli]/m0 - uSqr, 0.0);

    Meq_(0) = moments(0)[celli];
    Meq_(1) = moments(1)[celli];
    Meq_(2) = moments(2)[celli];
    Meq_(3) = m0*(3.0*sigma11*u + u*uSqr);
    Meq_(4) = m0*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);

    if (nDimensions_ > 1)
    {
        scalar v = meanVelocity(m0, moments(0,1)[celli]);
        scalar vSqr = sqr(v);
        scalar sigma22 = max(moments(0,2)[celli]/m0 - vSqr, 0.0);

        Meq_(0,1) = moments(0,1)[celli];
        Meq_(1,1) = moments(1,1)[celli];
        Meq_(0,2) = moments(0,2)[celli];
        Meq_(0,3) = m0*(3.0*sigma22*v + v*vSqr);
        Meq_(0,4) = m0*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
    }
    if (nDimensions_ > 2)
    {
        scalar w = meanVelocity(m0, moments(0,0,1)[celli]);
        scalar wSqr = sqr(w);
        scalar sigma33 = max(moments(0,0,2)[celli]/m0 - wSqr, 0.0);

        Meq_(0,0,1) = moments(0,0,1)[celli];
        Meq_(1,0,1) = moments(1,0,1)[celli];
        Meq_(0,1,1) = moments(0,1,1)[celli];
        Meq_(0,0,2) = moments(0,0,2)[celli];
        Meq_(0,0,3) = m0*(3.0*sigma33*w + w*wSqr);
        Meq_(0,0,4) = m0*(6.0*wSqr*sigma33 + 3.0*sqr(sigma33) + wSqr*wSqr);
    }
}

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateFields()
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    volScalarField m0(max(moments(0), SMALL));

    dimensionedScalar zeroVar("zero", sqr(dimVelocity), 0.0);

    volScalarField u(meanVelocity(m0, moments(1)));
    volScalarField uSqr(sqr(u));
    volScalarField sigma11(max(moments(2)/m0 - uSqr, zeroVar));
    Meqf_(0) = moments(0);
    Meqf_(1) = moments(1);
    Meqf_(2) = moments(2);
    Meqf_(3) = m0*(3.0*sigma11*u + u*uSqr);
    Meqf_(4) = m0*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);

    if (nDimensions_ > 1)
    {
        volScalarField v(meanVelocity(m0, moments(0,1)));
        volScalarField vSqr(sqr(v));
        volScalarField sigma22(max(moments(0,2)/m0 - vSqr, zeroVar));

        Meqf_(0,1) = moments(0,1);
        Meqf_(1,1) = moments(1,1);
        Meqf_(0,2) = moments(0,2);
        Meqf_(0,3) = m0*(3.0*sigma22*v + v*vSqr);
        Meqf_(0,4) = m0*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
    }
    if (nDimensions_ > 2)
    {
        volScalarField w(meanVelocity(m0, moments(0,0,1)));
        volScalarField wSqr(sqr(w));
        volScalarField sigma33(max(moments(0,0,2)/m0 - wSqr, zeroVar));

        Meqf_(0,0,1) = moments(0,0,1);
        Meqf_(1,0,1) = moments(1,0,1);
        Meqf_(0,1,1) = moments(0,1,1);
        Meqf_(0,0,2) = moments(0,0,2);
        Meqf_(0,0,3) = m0*(3.0*sigma33*w + w*wSqr);
        Meqf_(0,0,4) = m0*(6.0*wSqr*sigma33 + 3.0*sqr(sigma33) + wSqr*wSqr);
    }
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
    tauCollisional_("tau", dimTime, dict),
    Meqf_(quadrature.moments().size(), momentOrders_),
    Meq_(quadrature.moments().size(), momentOrders_)
{
    if (!ode)
    {
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

Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::explicitCollisionSource(const label mi, const label celli) const
{
    return (Meq_[mi] - quadrature_.moments()[mi][celli])/tauCollisional_.value();
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
