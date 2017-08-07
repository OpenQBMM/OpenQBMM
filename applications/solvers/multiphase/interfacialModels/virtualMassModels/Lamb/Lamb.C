/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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

#include "Lamb.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(Lamb, 0);
    addToRunTimeSelectionTable
    (
        virtualMassModel,
        Lamb,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::Lamb::Lamb
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    virtualMassModel(dict, pair, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModels::Lamb::~Lamb()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::virtualMassModels::Lamb::Cvm
(
    const label nodei,
    const label nodej
) const
{
    volScalarField E(min(max(pair_.E(nodei, nodej), SMALL), 1 - SMALL));
    volScalarField rtOmEsq(sqrt(1 - sqr(E)));

    return (rtOmEsq - E*acos(E))/(E*acos(E) - sqr(E)*rtOmEsq);
}


Foam::tmp<Foam::volScalarField> Foam::virtualMassModels::Lamb::Cvm() const
{
    volScalarField E(min(max(pair_.E(), SMALL), 1 - SMALL));
    volScalarField rtOmEsq(sqrt(1 - sqr(E)));

    return (rtOmEsq - E*acos(E))/(E*acos(E) - sqr(E)*rtOmEsq);
}


// ************************************************************************* //