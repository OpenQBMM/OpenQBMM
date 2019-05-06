/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 Alberto Passalacqua
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
            moments_(0).dimensions(),
            0.0
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
            symmTensor(-1.0)
        )
    ),
    isotropic_(true),
    nVelocityDims_(0)
{
    if (Sigma_[0].xx() > 0)
    {
        isotropic_ = false;
    }
    else
    {
        Sigma_ = symmTensor::I*Theta_;
    }
    labelList zeroOrder(momentOrders[0].size(), 0);
    forAll(momentOrders[0], cmpt)
    {
        labelList orderOne(zeroOrder);
        orderOne[cmpt] = 1;
        if
        (
            moments_(orderOne).dimensions()/moments_(0).dimensions()
         == dimVelocity
        )
        {
            nVelocityDims_++;
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
    const dictionary& dict
)
{
    if (dict.found("m0"))
    {
        m0_ = dimensionedScalar("m0", moments_(0).dimensions(), dict);
    }
    if (dict.found("U"))
    {
        U_ = dimensionedVector("U", dimVelocity, dict);
    }
    if (dict.found("Theta"))
    {
        Theta_ = dimensionedScalar("Theta", sqr(dimVelocity), dict);
    }
    if (dict.found("Sigma"))
    {
        Sigma_ = dimensionedSymmTensor("Sigma", sqr(dimVelocity), dict);
    }
}

void Foam::momentGenerationSubModels::gaussian::updateMoments
(
    const label celli
)
{
//     reset();
    if (isotropic_)
    {
        forAll (abscissae_[0], cmpt)
        {
            abscissae_[0][cmpt] = U_[celli][cmpt] + sqrt(Theta_[celli]);
            abscissae_[1][cmpt] = U_[celli][cmpt] - sqrt(Theta_[celli]);
        }
        weights_[0] = m0_[celli]/2.0;
        weights_[1] = m0_[celli]/2.0;
    }
    else
    {
        label nWeights = 2;
        abscissae_[0][0] = U_[celli].x() + sqrt(Sigma_[celli].xx());
        abscissae_[1][0] = U_[celli].x() - sqrt(Sigma_[celli].xx());

        if (nVelocityDims_ > 1)
        {
            nWeights = 4;

            abscissae_[0][1] = U_[celli].y() + sqrt(Sigma_[celli].yy());
            abscissae_[1][1] = U_[celli].y() - sqrt(Sigma_[celli].yy());

            abscissae_[2][0] = U_[celli].x() + sqrt(Sigma_[celli].xy());
            abscissae_[2][1] = U_[celli].y() + sqrt(Sigma_[celli].xy());
            abscissae_[3][0] = U_[celli].x() - sqrt(Sigma_[celli].xy());
            abscissae_[3][1] = U_[celli].y() - sqrt(Sigma_[celli].xy());
        }

        if (nVelocityDims_ > 2)
        {
            nWeights = 8;

            abscissae_[0][2] = U_[celli].z() + sqrt(Sigma_[celli].zz());
            abscissae_[1][2] = U_[celli].z() - sqrt(Sigma_[celli].zz());

            abscissae_[4][0] = U_[celli].x() + sqrt(Sigma_[celli].xz());
            abscissae_[4][2] = U_[celli].z() + sqrt(Sigma_[celli].xz());
            abscissae_[5][0] = U_[celli].x() - sqrt(Sigma_[celli].xz());
            abscissae_[5][2] = U_[celli].z() - sqrt(Sigma_[celli].xz());

            abscissae_[6][1] = U_[celli].y() + sqrt(Sigma_[celli].yz());
            abscissae_[6][2] = U_[celli].z() + sqrt(Sigma_[celli].yz());
            abscissae_[7][1] = U_[celli].y() - sqrt(Sigma_[celli].yz());
            abscissae_[7][2] = U_[celli].z() - sqrt(Sigma_[celli].yz());
        }

        for (label i = 0; i < nWeights; i++)
        {
            weights_[i] = m0_[celli]/scalar(nWeights);
        }
    }

    momentGenerationModel::updateMoments();
}

void Foam::momentGenerationSubModels::gaussian::updateMoments
(
    const label patchi,
    const label facei
)
{
    reset();
    const vector& u = U_.boundaryField()[patchi][facei];
    if (isotropic_)
    {
        forAll (abscissae_[0], cmpt)
        {
            abscissae_[0][cmpt] =
                u[cmpt] + sqrt(Theta_.boundaryField()[patchi][facei]);
            abscissae_[0][cmpt] =
                u[cmpt] - sqrt(Theta_.boundaryField()[patchi][facei]);
        }
        weights_[0] = m0_.boundaryField()[patchi][facei]/2.0;
        weights_[1] = m0_.boundaryField()[patchi][facei]/2.0;
    }
    else
    {
        label nWeights = 2;
        const symmTensor& sigma = Sigma_.boundaryField()[patchi][facei];
        abscissae_[0][0] = u.x() + sqrt(sigma.xx());
        abscissae_[1][0] = u.x() - sqrt(sigma.xx());

        if (nVelocityDims_ > 1)
        {
            nWeights = 4;

            abscissae_[0][1] = u.y() + sqrt(sigma.yy());
            abscissae_[1][1] = u.y() - sqrt(sigma.yy());

            abscissae_[2][0] = u.x() + sqrt(sigma.xy());
            abscissae_[2][1] = u.y() + sqrt(sigma.xy());
            abscissae_[3][0] = u.x() - sqrt(sigma.xy());
            abscissae_[3][1] = u.y() - sqrt(sigma.xy());
        }

        if (nVelocityDims_ > 2)
        {
            nWeights = 8;

            abscissae_[0][2] = u.z() + sqrt(sigma.zz());
            abscissae_[1][2] = u.z() - sqrt(sigma.zz());

            abscissae_[4][0] = u.x() + sqrt(sigma.xz());
            abscissae_[4][2] = u.z() + sqrt(sigma.xz());
            abscissae_[5][0] = u.x() - sqrt(sigma.xz());
            abscissae_[5][2] = u.z() - sqrt(sigma.xz());

            abscissae_[6][1] = u.y() + sqrt(sigma.yz());
            abscissae_[6][2] = u.z() + sqrt(sigma.yz());
            abscissae_[7][1] = u.y() - sqrt(sigma.yz());
            abscissae_[7][2] = u.z() - sqrt(sigma.yz());
        }

        for (label i = 0; i < nWeights; i++)
        {
            weights_[i] = m0_.boundaryField()[patchi][facei]/scalar(nWeights);
        }
    }
    momentGenerationModel::updateMoments();
}

// ************************************************************************* //
