/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 Alberto Passalacqua
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

    copyWeightsFromCells();

    const vectorField bfNorm(patch_.nf());

    const PtrList<volVelocityNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volVectorField& U = nodes[nodei].velocityAbscissae();

        vectorField& bfUOwn =
            nodesOwn_[nodei].velocityAbscissae().boundaryFieldRef()[patchi_];

        vectorField& bfUNei =
            nodesNei_[nodei].velocityAbscissae().boundaryFieldRef()[patchi_];

        const vectorField bfU
        (
            U.boundaryField()[patchi_].patchInternalField()
        );

        // Keep the abscissa of a node leaving the domain, and remove the
        // one of a node that would enter through it. The scaling once used
        // here was (bfU & bfSf), a volumetric flux, so the abscissa carried
        // the dimensions of a flux and scaled with the area of the face.
        bfUOwn = pos0(bfU & bfNorm)*bfU;
        bfUNei = bfUOwn;
    }
}


// ************************************************************************* //
