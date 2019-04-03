/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added additional return functions so that class can
                            be extended to polydisperse
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "vdfPhaseModel.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vdfPhaseModel::vdfPhaseModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    phaseModel(mesh, dict, phaseName),
    pbeDict_
    (
        IOobject
        (
            IOobject::groupName("populationBalanceProperties", phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    populationBalance_
    (
        populationBalanceModel::New
        (
            name_,
            pbeDict_,
            phiPtr_()
        )
    ),
    quadrature_
    (
        mesh.lookupObjectRef<velocityQuadratureApproximation>
        (
            IOobject::groupName
            (
                "quadratureProperties",
                phaseName
            )
        )
    ),
    computeVariance_(false),
    minD_
    (
        dimensionedScalar::lookupOrDefault
        (
            "minD",
            phaseDict_,
            dimLength,
            1e-6
        )
    )
{
    const labelList& velocityIndexes =
        quadrature_.nodes()[0].velocityIndexes();
    const labelListList& momentOrders =
        quadrature_.momentOrders();

    forAll(momentOrders, mi)
    {
        forAll(velocityIndexes, cmpt)
        {
            if (momentOrders[mi][velocityIndexes[cmpt]] >= 2)
            {
                computeVariance_ = true;
            }
        }
    }
    if (computeVariance_)
    {
         sigma_ = tmp<volTensorField>
         (
             new volTensorField
            (
                IOobject
                (
                    IOobject::groupName("Sigma", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedTensor("Sigma", sqr(dimVelocity), Zero)
            )
        );
        Theta_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Theta", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Theta", sqr(dimVelocity), 0.0)
            )
        );
    }

    correct();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vdfPhaseModel::~vdfPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::vdfPhaseModel::nNodes() const
{
    return quadrature_.nodes().size();
}

Foam::tmp<Foam::volScalarField> Foam::vdfPhaseModel::volumeFraction
(
    const label nodei
) const
{
    if (nodei == -1)
    {
        return quadrature_.moments()(0);
    }
    return quadrature_.nodes()[nodei].primaryWeight();
}

Foam::tmp<Foam::volScalarField>
Foam::vdfPhaseModel::d(const label nodei) const
{
    label sizeIndex = quadrature_.nodes()[0].sizeIndex();
    if (nodei == -1 || sizeIndex == -1)
    {
        return d_;
    }
    return
        Foam::max
        (
            quadrature_.nodes()[nodei].primaryAbscissae()[sizeIndex],
            minD_
        );

}

const Foam::volVectorField& Foam::vdfPhaseModel::U(const label nodei) const
{
    if (nodei == -1)
    {
        return U_;
    }
    return quadrature_.nodes()[nodei].velocityAbscissae();
}

Foam::volVectorField& Foam::vdfPhaseModel::U(const label nodei)
{
    if (nodei == -1)
    {
        return U_;
    }
    return quadrature_.nodes()[nodei].velocityAbscissae();
}

Foam::tmp<Foam::volVectorField>
Foam::vdfPhaseModel::Vs(const label nodei) const
{
    if (nodei == -1)
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("V", name_),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh(),
                dimensionedVector("zero", dimVelocity, Zero)
            )
        );
    }
    return U(nodei) - U_;
}

void Foam::vdfPhaseModel::solve()
{
    populationBalance_->solve();

    const labelList& velocityIndexes =
        quadrature_.nodes()[0].velocityIndexes();
    labelList orderZero(quadrature_.momentOrders()[0].size(), 0);
    volScalarField m0(quadrature_.moments()(0));
    volScalarField& alpha = *this;
    alpha = m0;
    m0.max(1e-10);

    forAll(velocityIndexes, cmpt)
    {
        labelList orderOne(orderZero);
        orderOne[velocityIndexes[cmpt]] = 1;

        volScalarField meanU(quadrature_.moments()(orderOne)/m0);
        U_.replace(cmpt, meanU);
    }
    phiPtr_() = fvc::flux(U_);
    alphaPhi_ = fvc::interpolate(*this)*phiPtr_();
    alphaRhoPhi_ = fvc::interpolate(rho())*alphaPhi_;

    label sizeIndex = quadrature_.nodes()[0].sizeIndex();
    if (sizeIndex != -1)
    {
        labelList orderOne(orderZero);
        orderOne[sizeIndex] = 1;
        d_ = quadrature_.moments()(orderOne)/Foam::max(alpha, residualAlpha_);
    }
}

void Foam::vdfPhaseModel::correct()
{
    quadrature_.updateMoments();

    const labelList& velocityIndexes =
        quadrature_.nodes()[0].velocityIndexes();
    labelList orderZero(quadrature_.momentOrders()[0].size(), 0);
    volScalarField m0(quadrature_.moments()(orderZero));
    m0.max(1e-10);

    forAll(velocityIndexes, cmpt)
    {
        labelList orderOne(orderZero);
        orderOne[velocityIndexes[cmpt]] = 1;

        volScalarField meanU(quadrature_.moments()(orderOne)/m0);
        U_.replace(cmpt, meanU);

        if (computeVariance_)
        {
            forAll(velocityIndexes, cmpt2)
            {
                labelList orderOne2(orderZero);
                labelList orderTwo(orderZero);
                orderOne2[velocityIndexes[cmpt2]] = 1;
                orderTwo[velocityIndexes[cmpt]] = 1;
                orderTwo[velocityIndexes[cmpt2]] += 1;

                volScalarField meanU2(quadrature_.moments()(orderOne2)/m0);

                volScalarField coVar
                (
                    quadrature_.moments()(orderTwo)/m0 - meanU*meanU2
                );
                sigma_.ref().replace(cmpt + cmpt2*3, coVar);
            }
        }
    }

    if (computeVariance_)
    {
        Theta_.ref() = tr(sigma_())/3.0;
    }

    label sizeIndex = quadrature_.nodes()[0].sizeIndex();
    volScalarField& alpha = *this;
    alpha = quadrature_.moments()(0);

    if (sizeIndex != -1)
    {
        labelList orderOne(quadrature_.momentOrders()[0].size(), 0);
        orderOne[sizeIndex] = 1;
        d_ = quadrature_.moments()(orderOne)/Foam::max(alpha, residualAlpha_);
    }
}


Foam::scalar Foam::vdfPhaseModel::realizableCo() const
{
   return populationBalance_->realizableCo();
}


Foam::scalar Foam::vdfPhaseModel::CoNum() const
{
    return populationBalance_->CoNum();
}

void Foam::vdfPhaseModel::relativeTransport()
{}


void Foam::vdfPhaseModel::averageTransport()
{}


void Foam::vdfPhaseModel::solveSource()
{}

// ************************************************************************* //
