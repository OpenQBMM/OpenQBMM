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

#include "reflectiveRotatingWallFvQuadraturePatch.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reflectiveRotatingWallFvQuadraturePatch, 0);

    addToRunTimeSelectionTable
    (
        fvQuadraturePatch,
        reflectiveRotatingWallFvQuadraturePatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflectiveRotatingWallFvQuadraturePatch
::reflectiveRotatingWallFvQuadraturePatch
(
    const fvPatch& patch,
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    PtrList<surfaceVelocityNode>& nodesOwn,
    PtrList<surfaceVelocityNode>& nodesNei
)
:
    reflectiveFvQuadraturePatch(patch, dict, quadrature, nodesOwn, nodesNei),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    omega_(Function1<scalar>::New("omega", dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reflectiveRotatingWallFvQuadraturePatch
::~reflectiveRotatingWallFvQuadraturePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::reflectiveRotatingWallFvQuadraturePatch
::wallTangentVelocity
(
    const vectorField& U,
    const vectorField& n
) const
{
    const fvMesh& mesh = quadrature_.nodes()[0].primaryWeight().mesh();
    const scalar t = mesh.time().timeOutputValue();
    scalar om = omega_->value(t);

    // Calculate the rotating wall velocity from the specification of the motion
    const vectorField Uw
    (
        (-om)*((patch_.Cf() - origin_) ^ (axis_/mag(axis_)))
    );

    // Remove the component of Up normal to the wall
    // just in case it is not exactly circular
    return (Uw - n*(n & Uw));
}


// ************************************************************************* //
