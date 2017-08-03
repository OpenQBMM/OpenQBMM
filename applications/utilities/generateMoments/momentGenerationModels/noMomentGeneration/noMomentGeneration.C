/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
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

#include "noMomentGeneration.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(noMomentGeneration, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        noMomentGeneration,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::noMomentGeneration::noMomentGeneration
(
    const dictionary& dict,
    const label nNodes,
    const bool extended,
    const bool Radau
)
:
    momentGenerationModel(dict, nNodes, extended, Radau)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::noMomentGeneration::~noMomentGeneration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::noMomentGeneration::updateQuadrature
(
    const dictionary& dict
)
{
    for (label mi = 0; mi < nMoments_; mi++)
    {
        moments_[mi] =
            dimensionedScalar
            (
                "moment."+ Foam::name(mi),
                moments_[mi].dimensions(),
                dict.lookup("moment."+ Foam::name(mi))
            );
    }
}


// ************************************************************************* //
