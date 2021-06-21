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

#include "alphaAndDiameterVelocity.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(alphaAndDiameterVelocity, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        alphaAndDiameterVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameterVelocity::alphaAndDiameterVelocity
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    alphaAndDiameter(mesh, dict, momentOrders, nNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameterVelocity::~alphaAndDiameterVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::alphaAndDiameterVelocity::updateMoments
(
    const dictionary& dict,
    const label patchi
)
{
    label size = reset();

    forAll(weights_, nodei)
    {
        const dictionary& nodeDict = dict.subDict("node" + Foam::name(nodei));
        vectorField U("U", nodeDict, size);

        for (label cmpt = 1; cmpt < abscissae_[nodei].size(); cmpt++)
        {
            abscissae_[nodei][cmpt] = U.component(cmpt - 1);
        }
    }

    alphaAndDiameter::updateMoments(dict, patchi);
}

void Foam::momentGenerationSubModels::alphaAndDiameterVelocity::updateMoments
(
    const dictionary& dict,
    const labelList& cells
)
{
    label size = reset();

    forAll(weights_, nodei)
    {
        const dictionary& nodeDict = dict.subDict("node" + Foam::name(nodei));
        vectorField U("U", nodeDict, size);

        for (label cmpt = 1; cmpt < abscissae_[nodei].size(); cmpt++)
        {
            abscissae_[nodei][cmpt] = U.component(cmpt - 1);
        }
    }

    alphaAndDiameter::updateMoments(dict, cells);
}

// ************************************************************************* //
