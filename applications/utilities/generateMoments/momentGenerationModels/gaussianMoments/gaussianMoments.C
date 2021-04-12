/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 Alberto Passalacqua
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

#include "gaussianMoments.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(gaussian, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        gaussian,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::gaussian::gaussian
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    momentGenerationModel(mesh, dict, momentOrders, nNodes),
    m0_
    (
        IOobject
        (
            IOobject::groupName("alpha", IOobject::group(dict.name())),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "m0",
            momentDimensions_(0),
            Zero
        )
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", IOobject::group(dict.name())),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector::lookupOrDefault
        (
            "U",
            dict,
            dimVelocity,
            Zero
        )
    ),
    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", IOobject::group(dict.name())),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar::lookupOrDefault
        (
            "Theta",
            dict,
            sqr(dimVelocity),
            -1.0
        )
    ),
    Sigma_
    (
        IOobject
        (
            IOobject::groupName("Sigma", IOobject::group(dict.name())),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor::lookupOrDefault
        (
            "Sigma",
            dict,
            sqr(dimVelocity),
            symmTensor::uniform(-1.0) // Placeholder (invalid value)
        )
    ),
    isotropic_(true),
    nVelocityDimensions_(0)
{
    if (Sigma_[0].xx() > 0)
    {
        isotropic_ = false;
    }
    else
    {
        Sigma_ = symmTensor::I*Theta_;
    }

    labelList zeroOrder(momentOrders[0].size(), Zero);
    
    forAll(momentOrders[0], cmpt)
    {
        labelList orderOne(zeroOrder);
        orderOne[cmpt] = 1;
    
        if (momentDimensions_(orderOne)/momentDimensions_(0) == dimVelocity)
        {
            nVelocityDimensions_++;
        }
    }

    if (isotropic_ && Theta_[0] < 0)
    {
        FatalErrorInFunction
            << "Granular temperature must be specified if velocity" << nl
            << "distribution is isotropic, and must be >= 0."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::gaussian::~gaussian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::gaussian::setNodes
(
    const dictionary& dict,
    const scalarField& alpha,
    const vectorField& U,
    const scalarField& Theta,
    const symmTensorField& Sigma
)
{
    label nWeights = 2;

    if (isotropic_)
    {
        forAll (abscissae_[0], cmpt)
        {
            abscissae_[0][cmpt] =
                U.component(cmpt) + sqrt(Theta);
            abscissae_[1][cmpt] =
                U.component(cmpt) - sqrt(Theta);
        }
    }
    else
    {
        label xx = symmTensor::XX;

        abscissae_[0][0] = U.component(0) + sqrt(Sigma.component(xx));
        abscissae_[1][0] = U.component(0) - sqrt(Sigma.component(xx));

        if (nVelocityDimensions_ > 1)
        {
            nWeights = 4;

            label yy = symmTensor::YY;
            label xy = symmTensor::XY;

            abscissae_[0][1] =
                U.component(1) + sqrt(Sigma.component(yy));
            abscissae_[1][1] =
                U.component(1) - sqrt(Sigma.component(yy));

            abscissae_[2][0] =
                U.component(0) + sqrt(Sigma.component(xy));
            abscissae_[2][1] =
                U.component(1) + sqrt(Sigma.component(xy));
            abscissae_[3][0] =
                U.component(0) - sqrt(Sigma.component(xy));
            abscissae_[3][1] =
                U.component(1) - sqrt(Sigma.component(xy));
        }

        if (nVelocityDimensions_ > 2)
        {
            nWeights = 8;

            label zz = symmTensor::ZZ;
            label xz = symmTensor::XZ;
            label yz = symmTensor::YZ;

            abscissae_[0][2] = U.component(2) + sqrt(Sigma.component(zz));
            abscissae_[1][2] = U.component(2) - sqrt(Sigma.component(zz));

            abscissae_[4][0] = U.component(0) + sqrt(Sigma.component(xz));
            abscissae_[4][2] = U.component(2) + sqrt(Sigma.component(xz));
            abscissae_[5][0] = U.component(0) - sqrt(Sigma.component(xz));
            abscissae_[5][2] = U.component(2) - sqrt(Sigma.component(xz));

            abscissae_[6][1] = U.component(1) + sqrt(Sigma.component(yz));
            abscissae_[6][2] = U.component(2) + sqrt(Sigma.component(yz));
            abscissae_[7][1] = U.component(1) - sqrt(Sigma.component(yz));
            abscissae_[7][2] = U.component(2) - sqrt(Sigma.component(yz));
        }
    }

    for (label i = 0; i < nWeights; i++)
    {
        weights_[i] = alpha/scalar(nWeights);
    }

    momentGenerationModel::updateMoments();
}

void Foam::momentGenerationSubModels::gaussian::updateMoments
(
    const dictionary& dict,
    const label patchi
)
{
    label size = reset(patchi);

    scalarField alpha;
    vectorField U;
    scalarField Theta;
    symmTensorField Sigma;

    if (dict.found("m0"))
    {
        alpha = scalarField("m0", dict, size);
    }
    else if (patchi == -1)
    {
        alpha = m0_.primitiveField();
    }
    else
    {
        alpha = m0_.boundaryField()[patchi];
    }

    if (dict.found("U"))
    {
        U = vectorField("U", dict, size);
    }
    else if (patchi == -1)
    {
        U = U_.primitiveField();
    }
    else
    {
        U = U_.boundaryField()[patchi];
    }

    if (dict.found("Theta"))
    {
        Theta = scalarField("Theta", dict, size);
    }
    else if (patchi == -1)
    {
        Theta = Theta_.primitiveField();
    }
    else
    {
        Theta = Theta_.boundaryField()[patchi];
    }

    if (dict.found("Sigma"))
    {
        Sigma = symmTensorField("Sigma", dict, size);
    }
    else if (patchi == -1)
    {
        Sigma = Sigma_.primitiveField();
    }
    else
    {
        Sigma = Sigma_.boundaryField()[patchi];
    }

    setNodes(dict, alpha, U, Theta, Sigma);
}


void Foam::momentGenerationSubModels::gaussian::updateMoments
(
    const dictionary& dict,
    const labelList& cells
)
{
    label size = reset(cells);

    scalarField alpha(size, Zero);
    vectorField U(size, Zero);
    scalarField Theta(size, Zero);
    symmTensorField Sigma(size, Zero);

    if (dict.found("m0"))
    {
        alpha = scalarField("m0", dict, size);
    }
    else
    {
        forAll(cells, celli)
        {
            alpha[celli] = m0_[cells[celli]];
        }
    }

    if (dict.found("U"))
    {
        U = vectorField("U", dict, size);
    }
    else
    {
        forAll(cells, celli)
        {
            U[celli] = U_[cells[celli]];
        }
    }

    if (dict.found("Theta"))
    {
        Theta = scalarField("Theta", dict, size);
    }
    else
    {
        forAll(cells, celli)
        {
            Theta[celli] = Theta_[cells[celli]];
        }
    }

    if (dict.found("Sigma"))
    {
        Sigma = symmTensorField("Sigma", dict, size);
    }
    else
    {
        forAll(cells, celli)
        {
            Sigma[celli] = Sigma_[cells[celli]];
        }
    }

    setNodes(dict, alpha, U, Theta, Sigma);
}

// ************************************************************************* //
