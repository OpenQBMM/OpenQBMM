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

#include "weightsAndAbscissae.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(weightsAndAbscissae, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        weightsAndAbscissae,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::weightsAndAbscissae
::weightsAndAbscissae
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

Foam::momentGenerationSubModels::weightsAndAbscissae
::~weightsAndAbscissae()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::weightsAndAbscissae::updateQuadrature
(
    const dictionary& dict
)
{
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        if (dict.found("node" + Foam::name(nodei)))
        {
            dictionary nodeDict(dict.subDict("node" + Foam::name(nodei)));
            if (nodei == 0 && radau_)
            {
                abscissae_[nodei].value() = 0.0;
            }
            else
            {
                abscissae_[nodei] = nodeDict.lookup("abscissa");
            }

            weights_[nodei] = nodeDict.lookup("weight");
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
