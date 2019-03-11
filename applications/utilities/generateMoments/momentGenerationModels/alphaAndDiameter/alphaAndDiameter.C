/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
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

#include "alphaAndDiameter.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(alphaAndDiameter, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        alphaAndDiameter,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameter::alphaAndDiameter
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    momentGenerationModel(mesh, dict, momentOrders, nNodes),
    alpha_
    (
        IOobject
        (
            IOobject::groupName
            (
                "alpha",
                IOobject::group(dict.name())
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        1.0
    ),
    scale_(dict.lookupOrDefault("scale", true)),
    rho_
    (
        IOobject
        (
            IOobject::groupName
            (
                "rho",
                IOobject::group(dict.name())
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar::lookupOrDefault("rho", dict, dimDensity, 0.0)
    ),
    ds_(nNodes, 0.0),
    alphas_(nNodes, 0.0),
    sumAlpha_(0.0)
{
    if (!dict.found("rho"))
    {
        autoPtr<rhoThermo> thermo = rhoThermo::New(mesh, alpha_.group());
        rho_ = thermo->rho();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameter::~alphaAndDiameter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::alphaAndDiameter::setNodes
(
    const dictionary& dict
)
{
    sumAlpha_ = 0.0;
    forAll(weights_, nodei)
    {
        word nodeName = "node" + Foam::name(nodei);
        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            ds_[nodei] = nodeDict.lookupType<scalar>("dia");
            alphas_[nodei] = nodeDict.lookupType<scalar>("alpha");
            sumAlpha_ += alphas_[nodei];
        }
        else
        {
            ds_[nodei] = 0.0;
            alphas_[nodei] = 0.0;
        }
    }
    sumAlpha_ = max(sumAlpha_, 1e-8);
}

void Foam::momentGenerationSubModels::alphaAndDiameter::updateMoments
(
    const label celli
)
{
    reset();

    forAll(weights_, nodei)
    {
        scalar alpha = alpha_[celli]*alphas_[nodei];
        if (scale_)
        {
            alpha /= sumAlpha_;
        }

        scalar rho = rho_[celli];

        abscissae_[nodei][0] =
            Foam::constant::mathematical::pi/6.0*rho*pow3(ds_[nodei]);

        if (abscissae_[nodei][0] > SMALL)
        {
            weights_[nodei] = rho*alpha/abscissae_[nodei][0];
        }
        else
        {
            weights_[nodei] = 0.0;
        }
    }

    momentGenerationModel::updateMoments();
}

void Foam::momentGenerationSubModels::alphaAndDiameter::updateMoments
(
    const label patchi,
    const label facei
)
{
    reset();

    forAll(weights_, nodei)
    {
        scalar alpha =
            alpha_.boundaryField()[patchi][facei]*alphas_[nodei];
        if (scale_)
        {
            alpha /= sumAlpha_;
        }

        scalar rho = rho_.boundaryField()[patchi][facei];

        abscissae_[nodei][0] =
            Foam::constant::mathematical::pi/6.0*rho*pow3(ds_[nodei]);

        if (abscissae_[nodei][0] > SMALL)
        {
            weights_[nodei] = rho*alpha/abscissae_[nodei][0];
        }
        else
        {
            weights_[nodei] = 0.0;
        }
    }

    momentGenerationModel::updateMoments();
}

// ************************************************************************* //
