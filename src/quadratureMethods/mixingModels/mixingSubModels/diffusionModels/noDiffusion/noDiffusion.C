/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "noDiffusion.H"
#include "addToRunTimeSelectionTable.H"

#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace diffusionModels
{
    defineTypeNameAndDebug(noDiffusion, 0);

    addToRunTimeSelectionTable
    (
        diffusionModel,
        noDiffusion,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::diffusionModels::noDiffusion
::noDiffusion
(
    const dictionary& dict
)
:
    diffusionModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::diffusionModels::noDiffusion
::~noDiffusion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::mixingSubModels::diffusionModels::noDiffusion
::momentDiff
(
    const volScalarField& moment
) const
{
    tmp<volScalarField> noDiff
    (
        new volScalarField
        (
            IOobject
            (
                "noDiff",
                moment.mesh().time().timeName(),
                moment.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            moment.mesh(),
            dimensionedScalar("zero", inv(dimTime), 0.0)
        )
    );

    return fvm::Sp(noDiff, moment);
}

// ************************************************************************* //
