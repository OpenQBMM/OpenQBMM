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

        // The node leaving the domain keeps its abscissa, and the one
        // entering carries it mirrored about the plane of the patch
        bfUOwn = U.boundaryField()[patchi_].patchInternalField();
        bfUNei = bfUOwn - 2.0*(bfUOwn & bfNorm)*bfNorm;
    }
}


// ************************************************************************* //
