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


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::updateCells(const label celli)
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    scalar m0 = max(moments(0)[celli], SMALL);

    // Mean velocity
    scalar u = meanVelocity(m0, moments(1)[celli]);
    scalar v = 0.0;
    scalar w = 0.0;

    scalar uSqr = sqr(u);
    scalar vSqr = 0.0;
    scalar wSqr = 0.0;

    // Variances of velocities
    scalar sigma1 = max(moments(2)[celli]/m0 - uSqr, 0.0);
    scalar sigma2 = 0.0;
    scalar sigma3 = 0.0;
    Theta_[celli] = sigma1;

    if (nDimensions_ > 1)
    {
        v = meanVelocity(m0, moments(0,1)[celli]);
        vSqr = sqr(v);
        sigma2 = max(moments(0,2)[celli]/m0 - vSqr, 0.0);
        Theta_[celli] += sigma2;
    }
    if (nDimensions_ > 2)
    {
        w = meanVelocity(m0, moments(0,0,1)[celli]);
        wSqr = sqr(w);
        sigma3 = max(moments(0,0,2)[celli]/m0 - wSqr, 0.0);
        Theta_[celli] += sigma3;
    }
    Theta_[celli] /= nDimensions_;

    scalar sigma11 = a1_*Theta_[celli] + b1_*sigma1;
    Meq_(0) = moments(0)[celli];
    Meq_(1) = moments(1)[celli];
    Meq_(2) = m0*(sigma11 + uSqr);
    Meq_(3) = m0*(3.0*sigma11*u + u*uSqr);
    Meq_(4) = m0*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);

    if (nDimensions_ > 1)
    {
        scalar sigma22 = a1_*Theta_[celli] + b1_*sigma2;
        scalar sigma12 = b1_*(moments(1,1)[celli]/m0 - u*v);
        Meq_(0,1) = moments(0,1)[celli];
        Meq_(1,1) = m0*(sigma12 + u*v);
        Meq_(0,2) = m0*(sigma22 + vSqr);
        Meq_(0,3) = m0*(3.0*sigma22*v + v*vSqr);
        Meq_(0,4) = m0*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
    }

    if (nDimensions_ > 2)
    {
        scalar sigma33 = a1_*Theta_[celli] + b1_*sigma3;
        scalar sigma13 = b1_*(moments(1,0,1)[celli]/m0 - u*w);
        scalar sigma23 = b1_*(moments(0,1,1)[celli]/m0 - v*w);
        Meq_(0,0,1) = moments(0,0,1)[celli];
        Meq_(1,0,1) = m0*(sigma13 + u*w);
        Meq_(0,1,1) = m0*(sigma23 + v*w);
        Meq_(0,0,2) = m0*(sigma33 + wSqr);
        Meq_(0,0,3) = m0*(3.0*sigma33*w + w*wSqr);
        Meq_(0,0,4) = m0*(6.0*wSqr*sigma33 + 3.0*sqr(sigma33) + wSqr*wSqr);
    }
}


void Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::updateFields()
{
    const volVectorMomentFieldSet& moments = quadrature_.moments();
    volScalarField m0(max(moments(0), SMALL));

    // Mean velocity
    volScalarField u(meanVelocity(m0, moments(1)));
    tmp<volScalarField> vTmp;
    tmp<volScalarField> wTmp;

    volScalarField uSqr(sqr(u));
    tmp<volScalarField> vSqrTmp;
    tmp<volScalarField> wSqrTmp;

    // Variances of velocities
    dimensionedScalar zeroVar("zero", sqr(dimVelocity), 0.0);
    volScalarField sigma1(max(moments(2)/m0 - uSqr, zeroVar));
    tmp<volScalarField> sigma2;
    tmp<volScalarField> sigma3;
    Theta_ = sigma1;

    if (nDimensions_ > 1)
    {
        vTmp = meanVelocity(m0, moments(0,1));
        vSqrTmp = sqr(vTmp());
        sigma2 = max(moments(0,2)/m0 - vSqrTmp(), zeroVar);
        Theta_ += sigma2();
    }
    if (nDimensions_ > 2)
    {
        wTmp = meanVelocity(m0, moments(0,0,1));
        wSqrTmp = sqr(wTmp());
        sigma3 = max(moments(0,0,2)/m0 - wSqrTmp, zeroVar);
        Theta_ += sigma3();
    }
    Theta_ /= nDimensions_;

    volScalarField sigma11(a1_*Theta_ + b1_*sigma1);

    Meqf_(0) = moments(0);
    Meqf_(1) = moments(1);
    Meqf_(2) = m0*(sigma11 + uSqr);
    Meqf_(3) = m0*(3.0*sigma11*u + u*uSqr);
    Meqf_(4) = m0*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);

    if (nDimensions_ > 1)
    {
        const volScalarField& v = vTmp();
        const volScalarField& vSqr = vSqrTmp();
        volScalarField sigma12(b1_*(moments(1,1)/m0 - u*v));
        volScalarField sigma22(a1_*Theta_ + b1_*sigma2);

        Meqf_(0,1) = moments(0,1);
        Meqf_(1,1) = m0*(sigma12 + u*v);
        Meqf_(0,2) = m0*(sigma22 + vSqr);
        Meqf_(0,3) = m0*(3.0*sigma22*v + v*vSqr);
        Meqf_(0,4) = m0*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);

        if (nDimensions_ > 2)
        {
            const volScalarField& w = wTmp();
            const volScalarField& wSqr = wSqrTmp();
            volScalarField sigma13(b1_*(moments(1,0,1)/m0 - u*w));
            volScalarField sigma23(b1_*(moments(0,1,1)/m0 - v*w));
            volScalarField sigma33(a1_*Theta_ + b1_*sigma3);

            Meqf_(0,0,1) = moments(0,0,1);
            Meqf_(1,0,1) = m0*(sigma13 + u*w);
            Meqf_(0,1,1) = m0*(sigma23 + v*w);
            Meqf_(0,0,2) = m0*(sigma33 + wSqr);
            Meqf_(0,0,3) = m0*(3.0*sigma33*w + w*wSqr);
            Meqf_(0,0,4) = m0*(6.0*wSqr*sigma33 + 3.0*sqr(sigma33) + wSqr*wSqr);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::esBGKCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature,
    const bool ode
)
:
    collisionKernel(dict, mesh, quadrature, ode),
    e_(dict.lookupType<scalar>("e")),
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
        dimensionedScalar("0", sqr(dimVelocity), 0.0)
    ),
    zeta_(dict_.lookupOrDefault("zeta", 1.0)),
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
    Meqf_(quadrature.moments().size(), momentOrders_),
    Meq_(quadrature.moments().size(), momentOrders_)
{
    scalar omega = (1.0 + e_)/2.0;
    scalar gamma = 1.0 - b_;
    a1_ = gamma*sqr(omega);
    b1_ = a1_ - 2.0*gamma*omega + 1.0;

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

Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::~esBGKCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::explicitCollisionSource(const label mi, const label celli) const
{
    scalar c = quadrature_.moments()[0][celli]/0.63;
    scalar gs0 = (2.0 - c)/(2.0*pow3(1.0 - c)) + 1.1603*c;
    scalar tauC =
        zeta_*sqrt(Foam::constant::mathematical::pi)*dp_()[celli]
       /max
        (
            12.0*gs0*quadrature_.moments()[0][celli]*sqrt(Theta_[celli]),
            1e-10
        );

    return (Meq_[mi] - quadrature_.moments()[mi][celli])/tauC;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::collisionKernels::esBGKCollision
::implicitCollisionSource(const volVectorMoment& m) const
{
    volScalarField c(quadrature_.moments()[0]/0.63);
    volScalarField gs0((2.0 - c)/(2.0*pow3(1.0 - c)) + 1.1603*c);
    volScalarField tauC
    (
        zeta_*sqrt(Foam::constant::mathematical::pi)*dp_()
       /max
        (
            12.0*gs0*quadrature_.moments()[0]*sqrt(Theta_),
            dimensionedScalar("small", dimVelocity, 1e-10)
        )
    );

    return
    (
        Meqf_(m.cmptOrders())/tauC
      - fvm::Sp(1/tauC, m)
    );
}
// ************************************************************************* //
