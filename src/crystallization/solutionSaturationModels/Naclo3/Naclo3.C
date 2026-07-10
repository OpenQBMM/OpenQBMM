/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "Naclo3.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solutionSaturationModels
{
    defineTypeNameAndDebug(Naclo3, 0);
    addToRunTimeSelectionTable(solutionSaturationModel, Naclo3, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::Naclo3::Naclo3
(
    const dictionary& dict,
    const objectRegistry& db
)
:
    solutionSaturationModel(db)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::Naclo3::~Naclo3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::solutionSaturationModels::Naclo3::Csat
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tCsat
    (
        volScalarField::New
        (
            "Csat",
            T.mesh(),
            dimensionedScalar(dimless)
        )
    );

    volScalarField& Csat = tCsat.ref();

    forAll(Csat, celli)
    {
        Csat[celli] =  7.12889 + 0.08455*(T[celli] - 273.0);//；mol/kg
    }

    volScalarField::Boundary& CsatBf = Csat.boundaryFieldRef();

    forAll(Csat.boundaryField(), patchi)
    {
        scalarField& Csatp = CsatBf[patchi];
        const scalarField& pp = T.boundaryField()[patchi];

        forAll(Csatp, facei)
        {
            Csatp[facei] = 7.12889 + 0.08455*(pp[facei] - 273.0);
        }
    }

    return tCsat;
}


// ************************************************************************* //
