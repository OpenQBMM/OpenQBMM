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

#include "outflowFvQuadraturePatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(outflowFvQuadraturePatch, 0);

    addToRunTimeSelectionTable
    (
        fvQuadraturePatch,
        outflowFvQuadraturePatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outflowFvQuadraturePatch::outflowFvQuadraturePatch
(
    const fvPatch& patch,
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    PtrList<surfaceVelocityNode>& nodesOwn,
    PtrList<surfaceVelocityNode>& nodesNei
)
:
    fvQuadraturePatch(patch, dict, quadrature, nodesOwn, nodesNei)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::outflowFvQuadraturePatch::~outflowFvQuadraturePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outflowFvQuadraturePatch::update()
{
    if (!patch_.size())
    {
        return;
    }

    const PtrList<volVelocityNode>& nodes = quadrature_.nodes();
    const fvMesh& mesh = nodes[0].primaryWeight().mesh();

    const vectorField& bfSf(mesh.Sf().boundaryField()[patchi_]);
    vectorField bfNorm(bfSf/mag(bfSf));

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

        vectorField bfU(U.boundaryField()[patchi_].patchInternalField());
        vectorField Un(bfU/max(mag(bfU), SMALL));
        bfUOwn = Foam::max((bfU & bfSf), scalar(0))*Un;
        bfUNei = bfUOwn;
    }
}


// ************************************************************************* //
