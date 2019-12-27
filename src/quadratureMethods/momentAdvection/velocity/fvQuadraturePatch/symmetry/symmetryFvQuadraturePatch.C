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

#include "symmetryFvQuadraturePatch.H"
#include "symmetryFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symmetryFvQuadraturePatch, 0);

    addToRunTimeSelectionTable
    (
        fvQuadraturePatch,
        symmetryFvQuadraturePatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::symmetryFvQuadraturePatch::symmetryFvQuadraturePatch
(
    const fvPatch& patch,
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    PtrList<surfaceVelocityNode>& nodesOwn,
    PtrList<surfaceVelocityNode>& nodesNei
)
:
    fvQuadraturePatch(patch, dict, quadrature, nodesOwn, nodesNei)
{
    if (!isA<symmetryFvPatch>(patch_))
    {
        FatalErrorInFunction
            << "Symmetry physical boundary required, but "
            << patch_.type() << " specified."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::symmetryFvQuadraturePatch::~symmetryFvQuadraturePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::symmetryFvQuadraturePatch::update()
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

        bfUOwn = U.boundaryField()[patchi_].patchInternalField();
        bfUNei = (bfUOwn - 2.0*(bfUOwn & bfNorm)*bfNorm);
    }
}


// ************************************************************************* //
