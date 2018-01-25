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

#include "velocity.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(velocity, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        velocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::velocity::velocity
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    momentGenerationModel(dict, momentOrders, nNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::velocity::~velocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::velocity::updateQuadrature
(
    const dictionary& dict
)
{
    reset();
    forAll(weights_, nodei)
    {
        word nodeName = "node" + Foam::name(nodei);
        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            scalar dia(readScalar(nodeDict.lookup("dia")));
            scalar alpha(readScalar(nodeDict.lookup("alpha")));
            scalar rho(readScalar(nodeDict.lookup("rho")));
            vector U(nodeDict.lookup("U"));

            scalar mass =
                (4.0/3.0)*Foam::constant::mathematical::pi*rho*pow3(dia/2.0);

            weights_[nodei] = rho*alpha/mass;

            forAll(abscissae_[nodei], cmpti)
            {
                abscissae_[nodei][cmpti] = U.component(cmpti);
            }
        }
    }

    updateMoments();
}


// ************************************************************************* //
