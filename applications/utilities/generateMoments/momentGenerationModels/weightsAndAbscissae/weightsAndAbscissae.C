/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 Alberto Passalacqua
    Copyright (C) 2018 OpenFOAM Foundation
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
    const fvMesh& mesh,
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    momentGenerationModel(mesh, dict, momentOrders, nNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::weightsAndAbscissae
::~weightsAndAbscissae()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::weightsAndAbscissae::updateMoments
(
    const dictionary& dict,
    const label patchi
)
{
    label size = reset(patchi);
    tmp<scalarField> scale;

    if (dict.found("scale"))
    {
        scale = tmp<scalarField>(new scalarField("scale", dict, size));
    }

    forAll(weights_, nodei)
    {
        word nodeName = "node" + Foam::name(nodei);

        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));

            scalarList abscissae
            (
                nodeDict.found("abscissae")
              ? scalarList(nodeDict.lookup("abscissae"))
              : scalarList(momentOrders_[0].size(), 0.0)
            );

            forAll(abscissae, i)
            {
                forAll(abscissae_[nodei], cmpt)
                {
                    if (nodeDict.found("abscissae" + Foam::name(cmpt)))
                    {
                        abscissae_[nodei][cmpt] =
                            scalarField
                            (
                                "abscissae" + Foam::name(cmpt),
                                nodeDict,
                                size
                            );
                    }
                    else
                    {
                        abscissae_[nodei][cmpt] = abscissae[cmpt];
                    }
                }
            }

            weights_[nodei] = scalarField("weight", nodeDict, size);
            
            if (scale.valid())
            {
                weights_[nodei] *= scale();
            }
        }
    }

    // Update since weights and abscissae are constant
    momentGenerationModel::updateMoments();
}


void Foam::momentGenerationSubModels::weightsAndAbscissae::updateMoments
(
    const dictionary& dict,
    const labelList& cells
)
{
    label size = reset(cells);
    tmp<scalarField> scale;

    if (dict.found("scale"))
    {
        scale = tmp<scalarField>(new scalarField("scale", dict, size));
    }

    forAll(weights_, nodei)
    {
        word nodeName = "node" + Foam::name(nodei);

        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));

            scalarList abscissae
            (
                nodeDict.found("abscissae")
              ? scalarList(nodeDict.lookup("abscissae"))
              : scalarList(momentOrders_[0].size(), 0.0)
            );

            forAll(abscissae_[nodei], cmpt)
            {
                if (nodeDict.found("abscissae" + Foam::name(cmpt)))
                {
                    abscissae_[nodei][cmpt] =
                        scalarField
                        (
                            "abscissae" + Foam::name(cmpt),
                            nodeDict,
                            size
                        );
                }
                else
                {
                    abscissae_[nodei][cmpt] = abscissae[cmpt];
                }
            }

            weights_[nodei] = scalarField("weight", nodeDict, size);
            
            if (scale.valid())
            {
                weights_[nodei] *= scale();
            }
        }
    }

    // Update since weights and abscissae are constant
    momentGenerationModel::updateMoments();
}

// ************************************************************************* //
