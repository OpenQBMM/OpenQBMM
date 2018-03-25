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

#include "momentGenerationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentGenerationModel, 0);
    defineRunTimeSelectionTable(momentGenerationModel, dictionary);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationModel::updateMoments()
{
    for (label mi = 0; mi < nMoments_; mi++)
    {
        moments_[mi].value() = 0.0;

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            moments_[mi] += weights_[nodei]*pow(abscissae_[nodei], mi);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationModel::momentGenerationModel
(
    const dictionary& dict,
    const label nNodes,
    const bool extended,
    const bool radau
)
:
    dict_(dict),
    nNodes_((radau) ? (nNodes):(nNodes + 1)),
    nMoments_((extended) ? (2*nNodes + 1):(2*nNodes)),
    extended_(extended),
    radau_(radau),
    weights_(nNodes_),
    abscissae_(nNodes_),
    moments_(nMoments_)
{
    forAll(weights_, nodei)
    {
        weights_[nodei].dimensions().reset(dict.lookup("weightDimension"));
        abscissae_[nodei].dimensions().reset(dict.lookup("abscissaDimension"));
    }

    forAll(moments_, mi)
    {
        moments_[mi].dimensions().reset
        (
            (weights_[0]*pow(abscissae_[0], mi)).dimensions()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationModel::~momentGenerationModel()
{}


// ************************************************************************* //
