/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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

#include "Brownian.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(Brownian, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        Brownian,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::Brownian::Brownian
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    continuousPhase_(dict.lookupOrDefault("continuousPhase", word::null)),
    T_
    (
        dict.found("T")
      ? mesh.lookupObject<volScalarField>(dict.get<word>("T"))
      : mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("T", continuousPhase_)
        )
    ),
    mu_
    (
        dict.found("mu")
      ? mesh.lookupObject<volScalarField>(dict.get<word>("mu"))
      : mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("thermo:mu", continuousPhase_)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::Brownian::~Brownian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::populationBalanceSubModels::aggregationKernels::Brownian::Ka
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli,
    const label environment
) const
{
    return 2.0*Foam::constant::physicoChemical::k.value()*T_[celli]
            *sqr(d1 + d2)/(3.0*mu_[celli]
            *max(d1*d2, SMALL));
}

// ************************************************************************* //
