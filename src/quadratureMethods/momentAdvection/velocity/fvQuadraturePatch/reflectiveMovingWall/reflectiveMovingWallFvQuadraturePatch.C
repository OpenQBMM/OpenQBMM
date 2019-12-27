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

#include "reflectiveMovingWallFvQuadraturePatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reflectiveMovingWallFvQuadraturePatch, 0);

    addToRunTimeSelectionTable
    (
        fvQuadraturePatch,
        reflectiveMovingWallFvQuadraturePatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflectiveMovingWallFvQuadraturePatch
::reflectiveMovingWallFvQuadraturePatch
(
    const fvPatch& patch,
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    PtrList<surfaceVelocityNode>& nodesOwn,
    PtrList<surfaceVelocityNode>& nodesNei
)
:
    reflectiveFvQuadraturePatch(patch, dict, quadrature, nodesOwn, nodesNei),
    wallVelocity_("wallVelocity", dict, patch_.size())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reflectiveMovingWallFvQuadraturePatch
::~reflectiveMovingWallFvQuadraturePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::vectorField>
Foam::reflectiveMovingWallFvQuadraturePatch::wallTangentVelocity
(
    const vectorField& U,
    const vectorField& n
) const
{
    return wallVelocity_;
}

// ************************************************************************* //
