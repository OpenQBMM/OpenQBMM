/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
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

#include "esBGKCollision.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{
    defineTypeNameAndDebug(esBGKCollision, 0);

    addToRunTimeSelectionTable
    (
        collisionKernel,
        esBGKCollision,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * * * //

Foam::symmTensor
Foam::populationBalanceSubModels::collisionKernels::esBGKCollision::covariance
(
    const label celli,
    const scalar& u,
    const scalar& v,
    const scalar& w
)
{
    symmTensor sigma(Zero);

    const volVelocityMomentFieldSet& moments = quadrature_.moments();
    scalar m0 = max(moments(0)[celli], SMALL);

    // Variances of velocities
    scalar sigma1 = max(moments(2)[celli]/m0 - sqr(u), scalar(0));
    scalar sigma2 = 0.0;
    scalar sigma3 = 0.0;
    Theta_[celli] = sigma1;

    if (nDimensions_ > 1)
    {
        sigma2 = max(moments(0,2)[celli]/m0 - sqr(v), scalar(0));
        Theta_[celli] += sigma2;
    }

    if (nDimensions_ > 2)
    {
        sigma3 = max(moments(0,0,2)[celli]/m0 - sqr(w), scalar(0));
        Theta_[celli] += sigma3;
    }

    Theta_[celli] /= nDimensions_;
    sigma.xx() = a1_*Theta_[celli] + b1_*sigma1;

    if (nDimensions_ > 1)
    {
        sigma.yy() = a1_*Theta_[celli] + b1_*sigma2;
        sigma.xy() = b1_*(moments(1,1)[celli]/m0 - u*v);
    }

    if (nDimensions_ > 2)
    {
        sigma.zz() = a1_*Theta_[celli] + b1_*sigma3;
        sigma.xz() = b1_*(moments(1,0,1)[celli]/m0 - u*w);
        sigma.yz() = b1_*(moments(0,1,1)[celli]/m0 - v*w);
    }

    return sigma;
}


Foam::symmTensor
Foam::populationBalanceSubModels::collisionKernels::esBGKCollision::covariance
(
    const mappedScalarList& moments,
    const scalar& u,
    const scalar& v,
    const scalar& w
)
{
    symmTensor sigma(Zero);
    scalar m0 = max(moments(0), SMALL);

    // Variances of velocities
    scalar sigma1 = max(moments(2)/m0 - sqr(u), scalar(0));
    scalar sigma2 = 0.0;
    scalar sigma3 = 0.0;
    scalar Theta = sigma1;

    if (nDimensions_ > 1)
    {
        sigma2 = max(moments(0,2)/m0 - sqr(v), scalar(0));
        Theta += sigma2;
    }

    if (nDimensions_ > 2)
    {
        sigma3 = max(moments(0,0,2)/m0 - sqr(w), scalar(0));
        Theta += sigma3;
    }

    Theta /= nDimensions_;
    sigma.xx() = a1_*Theta + b1_*sigma1;

    if (nDimensions_ > 1)
    {
        sigma.yy() = a1_*Theta + b1_*sigma2;
        sigma.xy() = b1_*(moments(1,1)/m0 - u*v);
    }

    if (nDimensions_ > 2)
    {
        sigma.zz() = a1_*Theta + b1_*sigma3;
        sigma.xz() = b1_*(moments(1,0,1)/m0 - u*w);
        sigma.yz() = b1_*(moments(0,1,1)/m0 - v*w);
    }

    return sigma;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::esBGKCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature
)
:
    BGKCollision(dict, mesh, quadrature),
    e_(readScalar(dict.lookup("e"))),
    b_(dict.lookupOrDefault<scalar>("b", 0)),
    Theta_
    (
        IOobject
        (
            "esBGK:Theta",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", sqr(dimVelocity), Zero)
    ),
    zeta_(dict_.lookupOrDefault<scalar>("zeta", 1))
{
    scalar omega = (1.0 + e_)/2.0;
    scalar gamma = 1.0 - b_;
    a1_ = gamma*sqr(omega);
    b1_ = a1_ - 2.0*gamma*omega + 1.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::~esBGKCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
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

    if (nSizes_ > 0)
    {
        return Meq_(momentOrder)[celli];
    }

    scalar c = quadrature_.moments()(0)[celli]/0.63;
    scalar gs0 = (2.0 - c)/(2.0*pow3(1.0 - c)) + 1.1603*c;

    scalar tauC =
        zeta_*sqrt(Foam::constant::mathematical::pi)*dp_()[celli]
       /max
        (
            12.0*gs0*quadrature_.moments()[0][celli]*sqrt(Theta_[celli]),
            1e-10
        );

    return
        (
            Meq_(momentOrder)[celli]
          - quadrature_.moments()(momentOrder)[celli]
        )/tauC;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::implicitCollisionSource(const volVelocityMoment& m) const
{
    if (!implicit_)
    {
        return tmp<fvScalarMatrix>
        (
            new fvScalarMatrix
            (
                m,
                m.dimensions()*dimVolume/dimTime
            )
        );
    }

    volScalarField c(quadrature_.moments()[0]/0.63);
    volScalarField gs0((2.0 - c)/(2.0*pow3(1.0 - c)) + 1.1603*c);
    volScalarField tauC
    (
        zeta_*sqrt(Foam::constant::mathematical::pi)*dp_()
       /max
        (
            12.0*gs0*quadrature_.moments()[0]*sqrt(Theta_),
            dimensionedScalar("SMALL", dimVelocity, 1e-10)
        )
    );

    return
    (
        Meq_(m.cmptOrders())/tauC
      - fvm::Sp(1/tauC, m)
    );
}
// ************************************************************************* //
