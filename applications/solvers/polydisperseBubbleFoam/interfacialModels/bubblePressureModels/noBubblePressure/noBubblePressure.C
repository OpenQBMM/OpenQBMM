/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "noBubblePressure.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bubblePressureModels
{
    defineTypeNameAndDebug(noBubblePressure, 0);
    addToRunTimeSelectionTable
    (
        bubblePressureModel,
        noBubblePressure,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubblePressureModels::noBubblePressure::noBubblePressure
(
    const dictionary& dict,
    const phasePair& pair
)
:
    bubblePressureModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bubblePressureModels::noBubblePressure::~noBubblePressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::bubblePressureModels::noBubblePressure::Fi
(
    const label nodei,
    const label nodej
) const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "noBubblePressure:Fi",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedVector("zero", dimF, Zero)
        )
    );
}


Foam::tmp<Foam::volVectorField>
Foam::bubblePressureModels::noBubblePressure::F
(
    const label nodei,
    const label nodej
) const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "nobubblePressure:F",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedVector("zero", dimF, Zero)
        )
    );
}


// ************************************************************************* //
