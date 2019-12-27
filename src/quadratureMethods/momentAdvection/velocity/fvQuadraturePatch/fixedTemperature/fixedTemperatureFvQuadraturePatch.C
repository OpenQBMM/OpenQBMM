/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 Alberto Passalacqua
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

#include "fixedTemperatureFvQuadraturePatch.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fixedTemperatureFvQuadraturePatch, 0);

    addToRunTimeSelectionTable
    (
        fvQuadraturePatch,
        fixedTemperatureFvQuadraturePatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedTemperatureFvQuadraturePatch::fixedTemperatureFvQuadraturePatch
(
    const fvPatch& patch,
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    PtrList<surfaceVelocityNode>& nodesOwn,
    PtrList<surfaceVelocityNode>& nodesNei
)
:
    fvQuadraturePatch(patch, dict, quadrature, nodesOwn, nodesNei),
    wallTemperature_("T", dict, patch.size()),
    nVelocityCmpts_(quadrature_.nodes()[0].velocityIndexes().size()),
    order000_(quadrature.momentOrders()[0].size(), 0),
    order100_(order000_),
    order010_(order000_),
    order001_(order000_),
    order200_(order000_),
    order020_(order000_),
    order002_(order000_)
{
    if (!isA<wallFvPatch>(patch_))
    {
        FatalErrorInFunction
            << "Fixed temperature requires a wall type boundary, "
            << "but " << patch_.type() << " was specified."
            << abort(FatalError);
    }

    labelList velocityIndexes = quadrature.nodes()[0].velocityIndexes();

    order100_[velocityIndexes[0]] = 1;
    order200_[velocityIndexes[0]] = 2;

    if (nVelocityCmpts_ > 1)
    {
        order010_[velocityIndexes[1]] = 1;
        order020_[velocityIndexes[1]] = 2;
    }

    if (nVelocityCmpts_ > 2)
    {
        order001_[velocityIndexes[2]] = 1;
        order002_[velocityIndexes[2]] = 2;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fixedTemperatureFvQuadraturePatch::~fixedTemperatureFvQuadraturePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedTemperatureFvQuadraturePatch::update()
{
    if (!patch_.size())
    {
        return;
    }

    const fvMesh& mesh = nodesOwn_[0].primaryWeight().mesh();

    const volVelocityMomentFieldSet& moments = quadrature_.moments();

    const vectorField& bfSf(mesh.Sf().boundaryField()[patchi_]);
    vectorField bfNorm(bfSf/mag(bfSf));

    scalarField m0(max(moments(0).boundaryField()[patchi_], scalar(1e-8)));
    vectorField T(bfSf.size(), Zero);

    T.replace
    (
        0,
        max
        (
            moments(order200_).boundaryField()[patchi_]/m0
          - sqr(moments(order100_).boundaryField()[patchi_]/m0),
            scalar(1e-8)
        )
    );

    if (nVelocityCmpts_ > 1)
    {
        T.replace
        (
            1,
            max
            (
                moments(order020_).boundaryField()[patchi_]/m0
              - sqr(moments(order010_).boundaryField()[patchi_]/m0),
                scalar(1e-8)
            )
        );
    }

    if (nVelocityCmpts_ > 2)
    {
        T.replace
        (
            2,
            max
            (
                moments(order002_).boundaryField()[patchi_]/m0
              - sqr(moments(order001_).boundaryField()[patchi_]/m0),
                scalar(1e-8)
            )
        );
    }

    scalarField scale
    (
        sqrt(wallTemperature_*scalar(nVelocityCmpts_)
       /(T & vector(1.0, 1.0, 1.0)))
    );

    scalarField Gin(bfSf.size(), Zero);
    scalarField Gout(bfSf.size(), Zero);

    const PtrList<volVelocityNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volVelocityNode& node = nodes[nodei];
        surfaceVelocityNode& nodeNei(nodesNei_[nodei]);
        surfaceVelocityNode& nodeOwn(nodesOwn_[nodei]);

        const volScalarField& weight = node.primaryWeight();
        surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        surfaceScalarField& weightNei = nodeNei.primaryWeight();
        const volVectorField& U = node.velocityAbscissae();
        surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        surfaceVectorField& UNei = nodeNei.velocityAbscissae();

        scalarField& bfwOwn = weightOwn.boundaryFieldRef()[patchi_];
        scalarField& bfwNei = weightNei.boundaryFieldRef()[patchi_];
        vectorField& bfUOwn = UOwn.boundaryFieldRef()[patchi_];
        vectorField& bfUNei = UNei.boundaryFieldRef()[patchi_];

        bfwOwn = weight.boundaryField()[patchi_].patchInternalField();
        bfwNei = bfwOwn;

        bfUOwn = U.boundaryField()[patchi_].patchInternalField();
        bfUNei = (bfUOwn - 2.0*(bfUOwn & bfNorm)*bfNorm)*scale;

        Gin += max(scalar(0), bfUOwn & bfSf)*bfwOwn;
        Gout -= min(scalar(0), bfUNei & bfSf)*bfwNei;
    }

    scalarField weightScale(Gin/(Gout + SMALL));

    forAll(nodes, nodei)
    {
        nodesNei_[nodei].primaryWeight().boundaryFieldRef()[patchi_] *=
            weightScale;
    }
}


// ************************************************************************* //
