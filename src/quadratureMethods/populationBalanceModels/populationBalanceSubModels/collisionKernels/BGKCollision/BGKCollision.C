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
    if (nSizes_ <= 1)
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
    List<vector> Us;

    // Size conditioned covariance tensors
    List<symmTensor> Sigmas;

    // Volume fraction of sizes
    scalarList weights;

    label sizei = 0;

    // Construct conditional moments
    for (label cmi = 0; cmi < nSizes_; cmi++)
    {
        forAll(velocityMoments_, vmi)
        {
            velocityMoments_[sizei][vmi][celli] = 0.0;
        }

        labelList index(nodeIndexes_[0].size(), 0);
        index[sizeIndex_] = cmi;

        forAll(velocityNodeIndexes_, nodei)
        {
            forAll(velocityIndexes_, cmpt)
            {
                index[velocityIndexes_[cmpt]] =
                    velocityNodeIndexes_[nodei][cmpt];
            }

            const volVelocityNode& node = quadrature_.nodes()(index);

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

        scalar m0 = velocityMoments_[sizei](0)[celli];
        if (m0 > 1e-10)
        {
            weights.append(m0);

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
            Us.append(vector(u,v,w));
            Sigmas.append(covariance(velocityMoments_[sizei], celli, u, v, w));
            sizei++;
        }
    }

    label nSizes = Us.size();

    scalarListList Ks(nSizes, scalarList(nSizes, 0.0));
    scalar alphard = 0.0;
    for (label cmi = 0; cmi < nSizes; cmi++)
    {
        const volVelocityNode& node = quadrature_.nodes()(cmi);
        alphard +=
            weights[cmi]
           /max(node.primaryAbscissae()[sizeIndex_][celli], 1e-10);
    }

    scalar pi = Foam::constant::mathematical::pi;
    for (label cmi = 0; cmi < nSizes; cmi++)
    {
        const volVelocityNode& nodei = quadrature_.nodes()(cmi);

        scalar Thetai = max(tr(Sigmas[cmi])/3.0, 0.0);

        for (label cmj = 0; cmj < nSizes; cmj++)
        {
            const volVelocityNode& nodej = quadrature_.nodes()(cmj);

            scalar Thetaj = max(tr(Sigmas[cmj])/3.0, 0.0);

            scalar Thetaij = Thetai + Thetaj;
            symmTensor Sij =
                0.5*(Sigmas[cmi] + Sigmas[cmj] + symmTensor::one*Thetaij);

            scalar di = max(nodei.primaryAbscissae()[sizeIndex_][celli], 1e-10);
            scalar dj = max(nodej.primaryAbscissae()[sizeIndex_][celli], 1e-10);
            scalar dij = (di + dj)*0.5;
            scalar XiPow3 = pow3(dij/dj);
            scalar massi = pi/6.0*pow3(di)*rhop_()[celli];
            scalar massj = pi/6.0*pow3(dj)*rhop_()[celli];

            scalar muij = 2.0*massj/(massi + massj);
            scalar g0ij = 1.0/alphac + 3.0*di*dj*alphard/(sqr(alphac)*(di + dj));

            Ks[cmi][cmj] =
                24.0*g0ij*weights[cmj]*XiPow3*sqrt(Thetaij)/(sqrt(pi)*dij);

            symmTensor Sigmaij =
                Sigmas[cmi]
              + 0.5*(1.0 + e())*muij*(0.25*(1.0 + e())*muij*Sij - Sigmas[cmi]);

            vector Uij = Us[cmi];
            if (cmi != cmj)
            {
                Us[cmi] += 0.25*(1.0 + e())*muij*(Us[cmj] - Us[cmi]);
            }

            forAllIter(List<momentFunction>, equilibriumMomentFunctions_, iter)
            {
                (*iter)
                (
                    Gs_[cmi][cmj],
                    celli,
                    1.0,
                    Uij.x(),
                    Uij.y(),
                    Uij.z(),
                    Sigmaij
                );
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

        for (label cmi = 0; cmi < nSizes; cmi++)
        {
            for (label cmj = 0; cmj < nSizes; cmj++)
            {
                Meq_[mi][celli] +=
                    Ks[cmi][cmj]
                   *(
                        Gs_[cmi][cmj](order)[celli]*weights[cmi]
                      - velocityMoments_[cmi](order)[celli]
                    )
                   *pow
                    (
                        quadrature_.nodes()(cmi).primaryAbscissae()[sizeIndex_][celli],
                        si
                    );
            }
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
        dimensionedScalar::lookupOrDefault("tau", dict, dimTime, 0.0)
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
    Meq_(momentOrders_.size(), momentOrders_)
{
    if (nSizes_ > 1)
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

    if (nSizes_ > 1)
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
                                "conditinoedVelocityMoment" + Foam::name(cmi),
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
    if (nSizes_ > 1)
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
