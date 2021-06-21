/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 Alberto Passalacqua
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    scale_(dict.lookupOrDefault("scale", false)),
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
        dimensionedScalar("rho", dimDensity, 0.0)
    ),
    diameters_(nNodes),
    alphas_(nNodes),
    sumAlpha_(),
    massBased_(dict.lookupOrDefault("massBased", true))
{
    if (!dict.found("rho") && massBased_)
    {
        autoPtr<rhoThermo> thermo = rhoThermo::New(mesh, alpha_.group());
        rho_ = thermo->rho();
    }
    else
    {
        rho_.primitiveFieldRef() = scalarField("rho", dict, mesh.nCells());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameter::~alphaAndDiameter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::alphaAndDiameter::updateMoments
(
    const dictionary& dict,
    const label patchi
)
{
    label size = reset(patchi);
    sumAlpha_ = scalarField(size, Zero);
    forAll(weights_, nodei)
    {
        word nodeName = "node" + Foam::name(nodei);

        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            diameters_[nodei] = scalarField("dia", nodeDict, size);
            alphas_[nodei] = scalarField("alpha", nodeDict, size);
            sumAlpha_ += alphas_[nodei];
        }
        else
        {
            diameters_[nodei] = 0.0;
            alphas_[nodei] = 0.0;
        }
    }

    sumAlpha_ = max(sumAlpha_, SMALL);
    
    scalarField alpha
    (
        patchi == -1
      ? alpha_.primitiveField()
      : alpha_.boundaryField()[patchi]
    );
    
    scalarField rho
    (
        patchi == -1
      ? rho_.primitiveField()
      : rho_.boundaryField()[patchi]
    );

    forAll(weights_, nodei)
    {
        scalarField alphai(alpha*alphas_[nodei]);
        if (scale_)
        {
            alpha /= sumAlpha_;
        }

        if (massBased_)
        {
            abscissae_[nodei][0] =
                Foam::constant::mathematical::pi/6.0*rho*pow3(diameters_[nodei]);

            weights_[nodei] =
                pos(abscissae_[nodei][0] - SMALL)
               *alphai*rho/max(abscissae_[nodei][0], SMALL);
        }
        else
        {
            abscissae_[nodei][0] = diameters_[nodei];
            scalarField V(pow3(diameters_[nodei]));
            weights_[nodei] = pos(V - SMALL)*alphai/max(V, SMALL);
        }
    }

    momentGenerationModel::updateMoments();
}


void Foam::momentGenerationSubModels::alphaAndDiameter::updateMoments
(
    const dictionary& dict,
    const labelList& cells
)
{
    label size = reset(cells);
    sumAlpha_ = scalarField(size, Zero);
    
    forAll(weights_, nodei)
    {
        word nodeName = "node" + Foam::name(nodei);
        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            diameters_[nodei] = scalarField("dia", nodeDict, size);
            alphas_[nodei] = scalarField("alpha", nodeDict, size);
            sumAlpha_ += alphas_[nodei];
        }
        else
        {
            diameters_[nodei] = 0.0;
            alphas_[nodei] = 0.0;
        }
    }

    sumAlpha_ = max(sumAlpha_, SMALL);
    scalarField alpha(size, Zero);
    scalarField rho(size, Zero);
    
    forAll(cells, celli)
    {
        alpha[celli] = alpha_[cells[celli]];
        rho[celli] = rho_[cells[celli]];
    }

    forAll(weights_, nodei)
    {
        scalarField alphai(alpha*alphas_[nodei]);
        if (scale_)
        {
            alpha /= sumAlpha_;
        }

        if (massBased_)
        {
            abscissae_[nodei][0] =
                Foam::constant::mathematical::pi/6.0*rho*pow3(diameters_[nodei]);

            weights_[nodei] =
                pos(abscissae_[nodei][0] - SMALL)
               *alphai*rho/max(abscissae_[nodei][0], SMALL);
        }
        else
        {
            abscissae_[nodei][0] = diameters_[nodei];
            scalarField V(pow3(diameters_[nodei]));
            weights_[nodei] = pos(V - SMALL)*alphai/max(V, SMALL);
        }
    }

    momentGenerationModel::updateMoments();
}


// ************************************************************************* //
