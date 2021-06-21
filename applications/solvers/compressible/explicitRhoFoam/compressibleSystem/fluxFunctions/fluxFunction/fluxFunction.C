/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "fluxFunction.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"
#include "subCycle.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxFunction, 0);
    defineRunTimeSelectionTable(fluxFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFunction::fluxFunction(const fvMesh& mesh)
:
    mesh_(mesh),
    dict_(mesh.schemesDict().subDict("compressible")),
    own_
    (
        IOobject
        (
            "fluxFunction:own",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("own", dimless, 1.0)
    ),
    nei_
    (
        IOobject
        (
            "fluxFunction:nei",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("nei", dimless, -1.0)
    ),
    residualAlpha_
    (
        dimensionedScalar::lookupOrDefault
        (
            "alpha",
            dict_,
            dimless,
            0.0
        )
    ),
    cutoffMa_
    (
        dimensionedScalar::lookupOrDefault
        (
            "cutoffMa",
            dict_,
            dimVelocity,
            0.0
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxFunction::~fluxFunction()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
