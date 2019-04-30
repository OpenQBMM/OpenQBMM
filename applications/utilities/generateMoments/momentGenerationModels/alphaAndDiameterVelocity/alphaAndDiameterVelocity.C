/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Alberto Passalacqua
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
    alphaAndDiameter(mesh, dict, momentOrders, nNodes),
    Us_(nNodes, Zero)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameterVelocity::~alphaAndDiameterVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::alphaAndDiameterVelocity::setNodes
(
    const dictionary& dict
)
{
    alphaAndDiameter::setNodes(dict);
    forAll(weights_, nodei)
    {
        word nodeName = "node" + Foam::name(nodei);
        if(dict.found(nodeName))
        {
            Us_[nodei] = dict.subDict(nodeName).lookupType<vector>("U");
        }
        else
        {
            Us_[nodei] = Zero;
        }
    }
}

void Foam::momentGenerationSubModels::alphaAndDiameterVelocity::updateMoments
(
    const label celli
)
{
    reset();

    alphaAndDiameter::updateMoments(celli);
    forAll(weights_, nodei)
    {
        for (label cmpt = 1; cmpt < abscissae_[nodei].size(); cmpt++)
        {
            abscissae_[nodei][cmpt] = Us_[nodei][cmpt - 1];
        }
    }

    momentGenerationModel::updateMoments();
}

void Foam::momentGenerationSubModels::alphaAndDiameterVelocity::updateMoments
(
    const label patchi,
    const label facei
)
{
    reset();

    alphaAndDiameter::updateMoments(patchi, facei);
    forAll(weights_, nodei)
    {
        for (label cmpt = 1; cmpt < abscissae_[nodei].size(); cmpt++)
        {
            abscissae_[nodei][cmpt] = Us_[nodei][cmpt - 1];
        }
    }

    momentGenerationModel::updateMoments();
}

// ************************************************************************* //
