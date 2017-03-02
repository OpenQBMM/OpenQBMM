/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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

#include "univariateQuadratureMomentSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateQuadratureMomentSet::univariateQuadratureMomentSet
(
    const dictionary& dict,
    const label nMoments,
    const word& support,
    const label nFixedQuadraturePoints,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
:
    univariateMomentSet
    (
        nMoments,
        support,
        nFixedQuadraturePoints
    ),
    abscissae_(),
    weights_(),
    inverted_(false),
    minKnownAbscissa_(minKnownAbscissa),
    maxKnownAbscissa_(maxKnownAbscissa),
    momentInverter_
    (
        Foam::univariateMomentInversion::New(dict)
    )
{}

Foam::univariateQuadratureMomentSet::univariateQuadratureMomentSet
(
    const dictionary& dict,
    const scalarList& m,
    const word& support,
    const label nFixedQuadraturePoints,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
:
    univariateMomentSet
    (
        m,
        support,
        nFixedQuadraturePoints
    ),
    abscissae_(),
    weights_(),
    inverted_(false),
    minKnownAbscissa_(minKnownAbscissa),
    maxKnownAbscissa_(maxKnownAbscissa),
    momentInverter_
    (
        Foam::univariateMomentInversion::New(dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateQuadratureMomentSet::~univariateQuadratureMomentSet()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateQuadratureMomentSet::invert()
{
    momentInverter_().invert
    (
        (*this), weights_, abscissae_, minKnownAbscissa_, maxKnownAbscissa_
    );

    inverted_ = true;
}

void Foam::univariateQuadratureMomentSet::update()
{
    // Recomputing all the moments (even if they originally were not realizable)
    // from quadrature (projection step).
    univariateMomentSet::update(weights_, abscissae_);

    realizabilityChecked_ = false;
    inverted_ = false;
}

// ************************************************************************* //
