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

#include "LuoSvendsen.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(LuoSvendsen, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        LuoSvendsen,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen
::LuoSvendsen
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    breakupKernel(dict, mesh),
    continuousPhase_(dict.lookupOrDefault("continuousPhase", word::null)),
    Cb_("Cb", dimless, dict),
    epsilonExp_(readScalar(dict.lookup("epsilonExp"))),
    nuExp_(readScalar(dict.lookup("nuExp"))),
    sizeExp_(readScalar(dict.lookup("sizeExp"))),
    flTurb_
    (
        mesh_.lookupObject<turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                continuousPhase_
            )
        )
    ),
    epsilon_(flTurb_.epsilon()),
        mu_
    (
        dict.found("mu")
      ? mesh.lookupObject<volScalarField>(dict.get<word>("mu"))
      : mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("thermo:mu", continuousPhase_)
        )
    ),
    rho_
    (
        dict.found("rho")
      ? mesh.lookupObject<volScalarField>(dict.get<word>("rho"))
      : mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("rho", continuousPhase_)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen::~LuoSvendsen()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen::Kb
(
    const scalar& abscissa,
    const label celli,
    const label environment
) const
{
    return Cb_.value()*pow(epsilon_[celli], epsilonExp_)
        *pow(mu_[celli]/rho_[celli], nuExp_)*pow(abscissa, sizeExp_);
}

// ************************************************************************* //
