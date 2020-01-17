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
    symmTensor sigma(Zero);

    const volVelocityMomentFieldSet& moments = quadrature_.moments();
    scalar m0 = max(moments(0)[celli], SMALL);
    sigma.xx() = max(moments(2)[celli]/m0 - sqr(u), scalar(0));

    if (nDimensions_ > 1)
    {
        sigma.xy() = moments(1,1)[celli]/m0 - u*v;
        sigma.yy() = max(moments(0,2)[celli]/m0 - sqr(v), scalar(0));

        if (nDimensions_ > 2)
        {
            sigma.xz() = moments(1,0,1)[celli]/m0 - u*w;
            sigma.yz() = moments(0,1,1)[celli]/m0 - v*w;
            sigma.zz() = max(moments(0,0,2)[celli]/m0 - sqr(w), scalar(0));
        }
    }

    return sigma;
}

Foam::symmTensor
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::covariance
(
    const mappedScalarList& moments,
    const scalar& u,
    const scalar& v,
    const scalar& w
)
{
    symmTensor sigma(Zero);
    scalar m0 = moments(0);

    if (m0 < SMALL)
    {
        return sigma;
    }

    sigma.xx() = max(moments(2)/m0 - sqr(u), scalar(0));

    if (nDimensions_ > 1)
    {
        sigma.xy() = moments(1,1)/m0 - u*v;
        sigma.yy() = max(moments(0,2)/m0 - sqr(v), scalar(0));

        if (nDimensions_ > 2)
        {
            sigma.xz() = moments(1,0,1)/m0 - u*w;
            sigma.yz() = moments(0,1,1)/m0 - v*w;
            sigma.zz() = max(moments(0,0,2)/m0 - sqr(w), scalar(0));
        }
    }

    return sigma;
}


// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::BGKCollision
::updateCells(const label celli)
{
    //- Monodisperse case
    const volVelocityMomentFieldSet& moments = quadrature_.moments();

    scalar m0 = moments(0)[celli];

    if (nSizes_ == 0)
    {
        scalar u = moments(1)[celli]/max(m0, SMALL);
        scalar v = 0;
        scalar w = 0;

        if (nDimensions_ > 1)
        {
            v = moments(0,1)[celli]/max(m0, SMALL);

            if (nDimensions_ > 2)
            {
                w = moments(0,0,1)[celli]/max(m0, SMALL);
            }
        }

        symmTensor sigma = covariance(celli, u, v, w);

        // Temporay moments
        mappedScalarList tmpMoments(momentOrders_.size(), momentOrders_);

        // Compute equilibrium distribuition moments given mean velocity
        // and covariance tensor
        forAllIter(List<momentFunction>, equilibriumMomentFunctions_, iter)
        {
            (*iter)(tmpMoments, m0, u, v, w, sigma);
        }

        forAll(moments, mi)
        {
            Meq_[mi][celli] = tmpMoments[mi];
        }

        return;
    }

    // Reset moment sources
    forAll(Meq_, mi)
    {
        Meq_[mi][celli] = Zero;
    }

    // Construct conditional moments
    forAll(velocityMoments_, sizei)
    {
        Ks_[sizei] = scalarList(nSizes_, Zero);
        forAll(velocityMomentOrders_, mi)
        {
            velocityMoments_[sizei][mi] = Zero;
            forAll(velocityMoments_, sizej)
            {
                Gs_[sizei][sizej][mi] = Zero;
            }
        }
    }

    // Continuous phase volume fraction used for radial distribution
    scalar alphac = 1.0 - m0;

    // Particle diameters
    scalarList ds(nSizes_);

    // Size conditioned mean velocities
    List<vector> Us(nSizes_, Zero);

    // Size conditioned covariance tensors
    List<symmTensor> Sigmas(nSizes_, Zero);

    // Volume fraction of sizes
    scalarList weights(nSizes_, Zero);

    // Compute size conditioned velocity moments
    forAll(quadrature_.nodes(), nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        const volVelocityNode& node = quadrature_.nodes()[nodei];
        label sizei = nodeIndex[sizeIndex_];

        forAll(velocityMomentOrders_, mi)
        {
            const labelList& vOrder = velocityMomentOrders_[mi];

            scalar mCmpt = node.primaryWeight()[celli];

            forAll(vOrder, cmpt)
            {
                mCmpt *=
                    pow(node.velocityAbscissae()[celli][cmpt], vOrder[cmpt]);
            }
            
            velocityMoments_[sizei](vOrder) += mCmpt;
        }
    }

    // Set size conditioned mean velocities and covariance tensors
    forAll(velocityMoments_, sizei)
    {
        scalar m0i = velocityMoments_[sizei](0);
        weights[sizei] = m0i;

        ds[sizei] =
            max
            (
                quadrature_.nodes()(sizei).primaryAbscissae()[sizeIndex_][celli],
                minD_
            );

        // Only compute variance and mean if m0 is not SMALL
        if (m0i > minM0_)
        {
            Us[sizei].x() = velocityMoments_[sizei](1)/m0i;

            if (nDimensions_ > 1)
            {
                Us[sizei].y() = velocityMoments_[sizei](0,1)/m0i;

                if (nDimensions_ > 2)
                {
                    Us[sizei].z() = velocityMoments_[sizei](0,0,1)/m0i;
                }
            }

            Sigmas[sizei] =
                covariance
                (
                    velocityMoments_[sizei],
                    Us[sizei].x(),
                    Us[sizei].y(),
                    Us[sizei].z()
                );
        }
    }

    scalar alphard(0);

    forAll(velocityMoments_, sizei)
    {
        alphard += weights[sizei]/ds[sizei];
    }

    // Compute source term coefficients and pair equilibrium distributions
    forAll(velocityMoments_, sizei)
    {
        scalar di = ds[sizei];
        scalar Vi = Foam::constant::mathematical::pi/6.0*pow3(di);
        scalar ni= weights[sizei]/Vi;

        forAll(velocityMoments_, sizej)
        {
            scalar m0ij = weights[sizei] + weights[sizej];

            //- Do not compute null moment sets
            if (m0ij > minM0_)
            {
                scalar dj = ds[sizej];
                scalar Vj = Foam::constant::mathematical::pi/6.0*pow3(dj);
                scalar nj = weights[sizej]/Vj;
                scalar dij = (di + dj)*0.5;

                vector Uij = Us[sizei];
                symmTensor Sigmaij = Sigmas[sizei];

                scalar xi = Foam::constant::mathematical::pi*(ni*sqr(di) 
                    + nj*sqr(dj))/6.0;

                scalar g0ij =
                    1.0/alphac
                  + 1.5*xi*di*dj/(sqr(alphac)*(dij))
                  + 0.5*sqr(xi)/pow3(alphac)*sqr(di*dj/dij);

                if (sizei == sizej)
                {
                    // Total granular temperature
                    scalar Thetai = max(tr(Sigmas[sizei]), scalar(0));
                    symmTensor Si = Sigmas[sizei] + symmTensor::I*Thetai;

                    // Covariance tensor
                    Sigmaij +=
                        0.5*(1.0 + e())*(0.25*(1.0 + e())*Si - Sigmas[sizei]);

                    // Collision time scale
                    Ks_[sizei][sizej] =
                        24.0*g0ij*weights[sizej]*sqrt(Thetai)
                       /(sqrt(Foam::constant::mathematical::pi)*di);
                }
                else
                {
                    const mappedScalarList& vmi = velocityMoments_[sizei];
                    const mappedScalarList& vmj = velocityMoments_[sizej];

                    // Total granular temperature
                    scalar Thetaij =
                        max
                        (
                            (vmi(2) + vmj(2))/m0ij
                          - sqr((vmi(1) + vmj(1))/m0ij),
                            scalar(0)
                        );

                    if (nDimensions_ > 1)
                    {
                        Thetaij +=
                            max
                            (
                                (vmi(0,2) + vmj(0,2))/m0ij
                              - sqr((vmi(0,1) + vmj(0,1))/m0ij),
                                scalar(0)
                            );
                    }

                    if (nDimensions_ > 2)
                    {
                        Thetaij +=
                            max
                            (
                                (vmi(0,0,2) + vmj(0,0,2))/m0ij
                              - sqr((vmi(0,0,1) + vmj(0,0,1))/m0ij),
                                scalar(0)
                            );
                    }

                    Thetaij /= nDimensions_;

                    symmTensor Sij =
                        0.5
                       *(
                            Sigmas[sizei]
                          + Sigmas[sizej]
                          + symmTensor::one*Thetaij
                        );

                    scalar XiPow3 = pow3(dij/dj);
                    scalar massi = Vi*rhos_[sizei];
                    scalar massj = Vj*rhos_[sizej];
                    scalar muij = 2.0*massj/(massi + massj);

                    // Collision time scale calculation
                    Ks_[sizei][sizej] =
                        24.0*g0ij*weights[sizej]*XiPow3*sqrt(Thetaij)
                       /(sqrt(Foam::constant::mathematical::pi)*dij);

                    // Covariance
                    Sigmaij +=
                        0.5*(1.0 + e())
                       *muij*(0.25*(1.0 + e())*muij*Sij - Sigmas[sizei]);

                    // Mean velocity
                    Uij += 0.25*(1.0 + e())*muij*(Us[sizej] - Us[sizei]);
                }

                // Compute equilibrium distribuition moments given mean velocity
                // and covariance tensor
                forAllIter
                (
                    List<momentFunction>, equilibriumMomentFunctions_,
                    iter
                )
                {
                    (*iter)
                    (
                        Gs_[sizei][sizej],
                        1.0,
                        Uij.x(),
                        Uij.y(),
                        Uij.z(),
                        Sigmaij
                    );
                }
            }
        }
    }

    // Compute source terms
    forAll(momentOrders_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        label si = momentOrder[sizeIndex_];
        labelList order(velocityIndexes_.size(), 0);

        forAll(velocityIndexes_, cmpt)
        {
            order[cmpt] = momentOrder[velocityIndexes_[cmpt]];
        }

        scalar& Meqi = Meq_(momentOrder)[celli];

        forAll(velocityMoments_, sizei)
        {
            scalar mCmpt = 0;

            forAll(velocityMoments_, sizej)
            {
                mCmpt +=
                    Ks_[sizei][sizej]
                   *(
                        Gs_[sizei][sizej](order)*weights[sizei]
                      - velocityMoments_[sizei](order)
                    );
            }

            if (si > 0)
            {
                mCmpt *=
                    pow
                    (
                        quadrature_.nodes()(sizei).primaryAbscissae()[sizeIndex_][celli],
                        si
                    );
            }

            Meqi += mCmpt;
        }
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
        dimensionedScalar::lookupOrDefault("tau", dict, dimTime, SMALL)
    ),
    Meq_(momentOrders_.size(), momentOrders_),
    Ks_(nSizes_, scalarList(nSizes_, Zero)),
    minM0_(dict.lookupOrDefault("minM0", SMALL))
{
    if (nSizes_ > 0)
    {
        implicit_ = false;
    }

    forAll(momentOrders_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];

        Meq_.set
        (
            momentOrder,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "Meq",
                        mappedList<scalar>::listToWord(momentOrder)
                    ),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    quadrature.moments()[mi].dimensions(),
                    Zero
                )
            )
        );
    }

    // Add required equilibrium moment functions
    mappedLabelList map(velocityMomentOrders_.size(), velocityMomentOrders_, 0);

    addMomentFunction1(map, equilibriumMomentFunctions_, 0)

    addMomentFunction3(map, equilibriumMomentFunctions_, 0,0,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 0,1)
    addMomentFunction1(map, equilibriumMomentFunctions_, 1)

    addMomentFunction3(map, equilibriumMomentFunctions_, 0,0,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,1,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 0,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 1,0,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 1,1)
    addMomentFunction1(map, equilibriumMomentFunctions_, 2)

    addMomentFunction3(map, equilibriumMomentFunctions_, 0,0,3)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,1,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,2,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 0,3)
    addMomentFunction3(map, equilibriumMomentFunctions_, 1,0,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 1,1,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 1,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 2,0,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 2,1)
    addMomentFunction1(map, equilibriumMomentFunctions_, 3)

    addMomentFunction3(map, equilibriumMomentFunctions_, 0,0,4)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,1,3)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,2,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,3,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 0,4)
    addMomentFunction3(map, equilibriumMomentFunctions_, 1,0,3)
    addMomentFunction2(map, equilibriumMomentFunctions_, 1,3)
    addMomentFunction3(map, equilibriumMomentFunctions_, 2,0,2)
    addMomentFunction2(map, equilibriumMomentFunctions_, 2,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 3,0,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 3,1)
    addMomentFunction1(map, equilibriumMomentFunctions_, 4)

    addMomentFunction3(map, equilibriumMomentFunctions_, 0,0,5)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,1,4)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,2,3)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,3,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,4,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 0,5)
    addMomentFunction3(map, equilibriumMomentFunctions_, 1,0,4)
    addMomentFunction2(map, equilibriumMomentFunctions_, 1,4)
    addMomentFunction3(map, equilibriumMomentFunctions_, 2,0,3)
    addMomentFunction2(map, equilibriumMomentFunctions_, 2,3)
    addMomentFunction3(map, equilibriumMomentFunctions_, 3,0,2)
    addMomentFunction2(map, equilibriumMomentFunctions_, 3,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 4,0,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 4,1)
    addMomentFunction1(map, equilibriumMomentFunctions_, 5)

    addMomentFunction3(map, equilibriumMomentFunctions_, 0,1,5)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,2,4)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,4,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,5,1)
    addMomentFunction3(map, equilibriumMomentFunctions_, 1,0,5)
    addMomentFunction2(map, equilibriumMomentFunctions_, 1,5)
    addMomentFunction3(map, equilibriumMomentFunctions_, 2,0,4)
    addMomentFunction2(map, equilibriumMomentFunctions_, 2,4)
    addMomentFunction3(map, equilibriumMomentFunctions_, 4,0,2)
    addMomentFunction2(map, equilibriumMomentFunctions_, 4,2)
    addMomentFunction3(map, equilibriumMomentFunctions_, 5,0,1)
    addMomentFunction2(map, equilibriumMomentFunctions_, 5,1)

    addMomentFunction3(map, equilibriumMomentFunctions_, 0,2,5)
    addMomentFunction3(map, equilibriumMomentFunctions_, 0,5,1)
    addMomentFunction3(map, equilibriumMomentFunctions_, 2,0,5)
    addMomentFunction2(map, equilibriumMomentFunctions_, 2,5)
    addMomentFunction3(map, equilibriumMomentFunctions_, 5,0,2)
    addMomentFunction2(map, equilibriumMomentFunctions_, 5,2)

    if (nSizes_ > 0)
    {
        // Equilibrium conditional distributions
        velocityMoments_.resize(nSizes_);
        Gs_.resize(nSizes_);

        forAll(Gs_, cmi)
        {
            Gs_.set(cmi, new PtrList<mappedScalarList>(nSizes_));

            forAll(Gs_[cmi], cmj)
            {
                Gs_[cmi].set
                (
                    cmj,
                    new mappedScalarList
                    (
                        velocityMomentOrders_.size(),
                        velocityMomentOrders_,
                        Zero
                    )
                );
            }

            velocityMoments_.set
            (
                cmi,
                new mappedScalarList
                (
                    velocityMomentOrders_.size(),
                    velocityMomentOrders_,
                    Zero
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

    //- Return calculated source
    if (nSizes_ > 0)
    {
        return Meq_(momentOrder)[celli];
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
          - fvm::Sp(1.0/tauCollisional_, m)
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
