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


// * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * * * //

Foam::symmTensor
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::covariance
(
    const label celli,
    const scalar& u,
    const scalar& v,
    const scalar& w
)
{
    symmTensor sigma;

    const volVelocityMomentFieldSet& moments = quadrature_.moments();
    scalar m0 = max(moments(0)[celli], SMALL);

    sigma.xx() = max(moments(2)[celli]/m0 - sqr(u), 0.0);

    if (nDimensions_ > 1)
    {
        sigma.yy() = max(moments(0,2)[celli]/m0 - sqr(v), 0.0);
        if (nDimensions_ > 2)
        {
            sigma.zz() = max(moments(0,0,2)[celli]/m0 - sqr(w), 0.0);
        }
    }

    return sigma;
}

// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateCells(const label celli)
{
    const volVelocityMomentFieldSet& moments = quadrature_.moments();
    scalar m0 = max(moments(0)[celli], SMALL);

    // Mean velocity
    scalar u = moments(1)[celli]/m0;
    scalar v = 0.0;
    scalar w = 0.0;
    if (nDimensions_ > 0)
    {
        v = moments(0,1)[celli]/m0;

        if (nDimensions_ > 1)
        {
            w = moments(0,0,1)[celli]/m0;
        }
    }
    symmTensor sigma = covariance(celli, u, v, w);

    forAllIter(List<mf>, equilibriumMomentFunctions_, iter)
    {
        (*iter)(Meq_, celli, m0, u, v, w, sigma);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BGKCollision::BGKCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature
)
:
    collisionKernel(dict, mesh, quadrature),
    tauCollisional_
    (
        dimensionedScalar::lookupOrDefault("tau", dict, dimTime, 0.0)
    ),
    Meq_(velocityMomentOrders_.size(), velocityMomentOrders_)
{
    forAll(Meq_, mi)
    {
        const labelList& momentOrder = velocityMomentOrders_[mi];
        Meq_.set
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

    if (Meq_.found(0))
    {
        equilibriumMomentFunctions_.append(&moment000);
    }
    if (Meq_.found(1))
    {
        equilibriumMomentFunctions_.append(&moment100);
    }
    if (Meq_.found(2))
    {
        equilibriumMomentFunctions_.append(&moment200);
    }
    if (Meq_.found(3))
    {
        equilibriumMomentFunctions_.append(&moment300);
    }
    if (Meq_.found(4))
    {
        equilibriumMomentFunctions_.append(&moment400);
    }
    if (Meq_.found(5))
    {
        equilibriumMomentFunctions_.append(&moment500);
    }
    if (Meq_.found(0,1))
    {
        equilibriumMomentFunctions_.append(&moment010);
    }
    if (Meq_.found(0,2))
    {
        equilibriumMomentFunctions_.append(&moment020);
    }
    if (Meq_.found(0,3))
    {
        equilibriumMomentFunctions_.append(&moment030);
    }
    if (Meq_.found(0,4))
    {
        equilibriumMomentFunctions_.append(&moment040);
    }
    if (Meq_.found(0,5))
    {
        equilibriumMomentFunctions_.append(&moment050);
    }
    if (Meq_.found(0,0,1))
    {
        equilibriumMomentFunctions_.append(&moment001);
    }
    if (Meq_.found(0,0,2))
    {
        equilibriumMomentFunctions_.append(&moment002);
    }
    if (Meq_.found(0,0,3))
    {
        equilibriumMomentFunctions_.append(&moment003);
    }
    if (Meq_.found(0,0,4))
    {
        equilibriumMomentFunctions_.append(&moment004);
    }
    if (Meq_.found(0,0,5))
    {
        equilibriumMomentFunctions_.append(&moment005);
    }
    if (Meq_.found(1,1,0))
    {
        equilibriumMomentFunctions_.append(&moment110);
    }
    if (Meq_.found(1,1,0))
    {
        equilibriumMomentFunctions_.append(&moment101);
    }
    if (Meq_.found(1,0,1))
    {
        equilibriumMomentFunctions_.append(&moment011);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::~BGKCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::explicitCollisionSource
(
    const labelList& momentOrder,
    const label celli
) const
{
    if (implicit_)
    {
        return 0.0;
    }

    return
    (
        Meq_(momentOrder)[celli]
      - quadrature_.moments()(momentOrder)[celli]
    )/tauCollisional_.value();
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::implicitCollisionSource(const volVelocityMoment& m) const
{
    if (implicit_)
    {
        return
        (
            Meq_(m.cmptOrders())/tauCollisional_
          - fvm::Sp(1/tauCollisional_, m)
        );
    }

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
