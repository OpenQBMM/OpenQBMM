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

#include "alphaAndDiameter.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(alphaAndDiameter, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        alphaAndDiameter,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameter::alphaAndDiameter
(
    const dictionary& dict,
    const label nNodes,
    const bool extended,
    const bool radau
)
:
    momentGenerationModel(dict, nNodes, extended, radau)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::alphaAndDiameter::~alphaAndDiameter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::alphaAndDiameter::updateQuadrature
(
    const dictionary& dict
)
{
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {

        if (dict.found("node" + Foam::name(nodei)))
        {
            dictionary nodeDict(dict.subDict("node"+Foam::name(nodei)));
            dimensionedScalar dia(nodeDict.lookup("dia"));
            dimensionedScalar alpha(nodeDict.lookup("alpha"));
            dimensionedScalar rho(nodeDict.lookup("rho"));

            abscissae_[nodei]
                = (4.0/3.0)*Foam::constant::mathematical::pi*rho*pow3(dia/2.0);

            weights_[nodei] = rho*alpha/abscissae_[nodei];

            if (nodei == 0 && radau_)
            {
                abscissae_[nodei].value() = 0.0;
            }
        }
        else
        {
            abscissae_[nodei].value() = 0.0;
            weights_[nodei].value() = 0.0;
        }
    }

    updateMoments();
}


// ************************************************************************* //
