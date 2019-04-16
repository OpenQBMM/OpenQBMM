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
    symmTensor sigma(Zero);

    const volVelocityMomentFieldSet& moments = quadrature_.moments();
    scalar m0 = max(moments(0)[celli], small);
    sigma.xx() = max(moments(2)[celli]/m0 - sqr(u), 0.0);

    if (nDimensions_ > 1)
    {
        sigma.xy() = moments(1,1)[celli]/m0 - u*v;
        sigma.yy() = max(moments(0,2)[celli]/m0 - sqr(v), 0.0);

        if (nDimensions_ > 2)
        {
            sigma.xz() = moments(1,0,1)[celli]/m0 - u*w;
            sigma.yz() = moments(0,1,1)[celli]/m0 - v*w;
            sigma.zz() = max(moments(0,0,2)[celli]/m0 - sqr(w), 0.0);
        }
    }
    return sigma;
}


Foam::symmTensor
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::covariance
(
    const mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& u,
    const scalar& v,
    const scalar& w
)
{
    symmTensor sigma(Zero);
    scalar m0 = max(moments(0)[celli], small);
    sigma.xx() = max(moments(2)[celli]/m0 - sqr(u), 0.0);

    if (nDimensions_ > 1)
    {
        sigma.xy() = moments(1,1)[celli]/m0 - u*v;
        sigma.yy() = max(moments(0,2)[celli]/m0 - sqr(v), 0.0);

        if (nDimensions_ > 2)
        {
            sigma.xz() = moments(1,0,1)[celli]/m0 - u*w;
            sigma.yz() = moments(0,1,1)[celli]/m0 - v*w;
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

    //- Monodisperse case
    if (nSizes_ < 0)
    {
        scalar m0 = moments(0)[celli];

        scalar u = moments(1)[celli]/max(m0, small);
        scalar v = 0.0;
        scalar w = 0.0;
        if (nDimensions_ > 1)
        {
            v = moments(0,1)[celli]/max(m0, small);

            if (nDimensions_ > 2)
            {
                w = moments(0,0,1)[celli]/max(m0, small);
            }
        }
        symmTensor sigma = covariance(celli, u, v, w);

        forAllIter(List<momentFunction>, equilibriumMomentFunctions_, iter)
        {
            (*iter)(Meq_, celli, m0, u, v, w, sigma);
        }
        return;
    }

    scalar alpha = quadrature_.moments()(0)[celli];
//     scalar c = min(alpha/0.63, 0.999);
    scalar alphac = 1.0 - alpha;
//     scalar g0 = (2.0 - c)/(2.0*pow3(1.0 - c));

    // Size conditioned mean velocities
    List<vector> Us(nSizes_, Zero);

    // Size conditioned covariance tensors
    List<symmTensor> Sigmas(nSizes_, Zero);

    // Volume fraction of sizes
    scalarList weights(nSizes_, 0.0);

    // Construct conditional moments
    forAll(velocityMoments_, sizei)
    {
        forAll(velocityMomentOrders_, mi)
        {
            velocityMoments_[sizei][mi][celli] = 0.0;
            forAll(velocityMoments_, sizej)
            {
                Gs_[sizei][sizej][mi][celli] = 0.0;
            }
        }
    }

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
            velocityMoments_[sizei](vOrder)[celli] += mCmpt;
        }
    }

    forAll(velocityMoments_, sizei)
    {
        scalar m0 = velocityMoments_[sizei](0)[celli];
        weights[sizei] = m0;
        if (m0 > small)
        {
            scalar u = velocityMoments_[sizei](1)[celli]/m0;
            scalar v = 0.0;
            scalar w = 0.0;
            if (nDimensions_ > 1)
            {
                v = velocityMoments_[sizei](0,1)[celli]/m0;

                if (nDimensions_ > 2)
                {
                    w = velocityMoments_[sizei](0,0,1)[celli]/m0;
                }
            }
            Us[sizei] = vector(u,v,w);
            Sigmas[sizei] =
                covariance(velocityMoments_[sizei], celli, u, v, w);
        }
    }


    scalarListList Ks(nSizes_, scalarList(nSizes_, 0.0));
    scalar alphard = 0.0;
    forAll(velocityMoments_, sizei)
    {
        const volVelocityNode& node = quadrature_.nodes()(sizei);
        alphard +=
            weights[sizei]
           /max(node.primaryAbscissae()[sizeIndex_][celli], small);
    }

    scalar pi = Foam::constant::mathematical::pi;
    forAll(velocityMoments_, sizei)
    {
        const volVelocityNode& nodei = quadrature_.nodes()(sizei);
        scalar Thetai = max(tr(Sigmas[sizei])/3.0, 0.0);
        scalar di = nodei.primaryAbscissae()[sizeIndex_][celli];
        symmTensor Sigmai = Sigmas[sizei] + Thetai*symmTensor::I;

        forAll(velocityMoments_, sizej)
        {
            const volVelocityNode& nodej = quadrature_.nodes()(sizej);
            scalar dj = nodej.primaryAbscissae()[sizeIndex_][celli];
            if
            (
                weights[sizei] > small && di > small
             && weights[sizej] > small && dj > small
            )
            {
                scalar Thetaj = max(tr(Sigmas[sizej])/3.0, 0.0);

                scalar Thetaij = Thetai + Thetaj;
                symmTensor Sigmaij =
                    0.5
                   *(
                        Sigmas[sizei]
                      + Sigmas[sizej]
                      + symmTensor::one*Thetaij
                    );

                scalar dij = (di + dj)*0.5;
                scalar XiPow3 = pow3(dij/dj);
                scalar massi = pi/6.0*pow3(di);
                scalar massj = pi/6.0*pow3(dj);

                scalar muij = 2.0*massj/(massi + massj);
                scalar g0ij =
                    1.0/alphac
                  + 3.0*di*dj*alphard/(sqr(alphac)*(di + dj));

                Ks[sizei][sizej] =
                    24.0*g0ij*weights[sizej]*XiPow3*sqrt(Thetaij)
                   /(sqrt(pi)*dij);

                symmTensor Sigma =
                    Sigmai
                  + 0.5*(1.0 + e())
                   *muij*(0.25*(1.0 + e())*muij*Sigmaij - Sigmai);

                vector Uij = Us[sizei];
                if (sizei != sizej)
                {
                    Uij +=
                        0.25*(1.0 + e())*muij*(Us[sizej] - Us[sizei]);
                }

                forAllIter
                (
                    List<momentFunction>, equilibriumMomentFunctions_,
                    iter
                )
                {
                    (*iter)
                    (
                        Gs_[sizei][sizej],
                        celli,
                        1.0,
                        Uij.x(),
                        Uij.y(),
                        Uij.z(),
                        Sigma
                    );
                }
            }
        }
    }

    forAll(momentOrders_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        label si = momentOrder[sizeIndex_];
        labelList order(velocityIndexes_.size(), 0);
        forAll(velocityIndexes_, cmpt)
        {
            order[cmpt] = momentOrder[velocityIndexes_[cmpt]];
        }
        Meq_[mi][celli] = 0.0;

        forAll(velocityMoments_, sizei)
        {
            scalar mCmpt = 0.0;
            forAll(velocityMoments_, sizej)
            {
                mCmpt +=
                    Ks[sizei][sizej]
                   *(
                        Gs_[sizei][sizej](order)[celli]*weights[sizei]
                      - velocityMoments_[sizei](order)[celli]
                    );
            }
            mCmpt *=
                pow
                (
                    quadrature_.nodes()(sizei).primaryAbscissae()[sizeIndex_][celli],
                    si
                );
            Meq_[mi][celli] += mCmpt;
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
        dimensionedScalar::lookupOrDefault("tau", dict, dimTime, small)
    ),
    Meq_(momentOrders_.size(), momentOrders_)
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
                    0.0
                )
            )
        );
    }
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
            Gs_.set(cmi, new PtrList<mappedPtrList<volScalarField>>(nSizes_));
            forAll(Gs_[cmi], cmj)
            {
                Gs_[cmi].set
                (
                    cmj,
                    new mappedPtrList<volScalarField>
                    (
                        velocityMomentOrders_.size(),
                        velocityMomentOrders_
                    )
                );

                forAll(velocityMomentOrders_, mi)
                {
                    const labelList& momentOrder = velocityMomentOrders_[mi];
                    Gs_[cmi][cmj].set
                    (
                        momentOrder,
                        new volScalarField
                        (
                            IOobject
                            (
                                IOobject::groupName
                                (
                                    "GaussianMoment"
                                  + Foam::name(cmi)
                                  + Foam::name(cmj),
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
                                0.0
                            )
                        )
                    );
                }
            }

            velocityMoments_.set
            (
                cmi,
                new mappedPtrList<volScalarField>
                (
                    velocityMomentOrders_.size(),
                    velocityMomentOrders_
                )
            );
            forAll(velocityMomentOrders_, mi)
            {
                const labelList& momentOrder = velocityMomentOrders_[mi];
                velocityMoments_[cmi].set
                (
                    momentOrder,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "conditionedVelocityMoment" + Foam::name(cmi),
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
                            0.0
                        )
                    )
                );
            }
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
