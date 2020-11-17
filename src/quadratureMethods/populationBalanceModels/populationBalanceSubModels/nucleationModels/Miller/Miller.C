/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
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

#include "Miller.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(Miller, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        Miller,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::Miller::Miller
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    continuousPhase_(dict.lookupOrDefault("continuousPhase", word::null)),
    MCarbon_("MCarbon", dimMass/dimMoles, dict),
    nCarbonDimer_("nCarbonDimer", dimless, dict),
    nCarbonPAM_("nCarbonPAM", dimless, dict),
    rhoSoot_("rhoSoot", dimDensity, dict),
    pamConcentration_
    (
        IOobject
        (
            "YPAM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_
    (
        dict.found("T")
      ? mesh.lookupObject<volScalarField>(dict.get<word>("T"))
      : mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("T", continuousPhase_)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::Miller::~Miller()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::Miller::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    scalar NA = Foam::constant::physicoChemical::NA.value();
    scalar MCarbon = MCarbon_.value();

    scalar abscissaNucleation =
        2.0*MCarbon*nCarbonDimer_.value()/(rhoSoot_.value()*NA);

    return 4.4*sqrt(Foam::constant::mathematical::pi
        *Foam::constant::physicoChemical::k.value()*T_[celli]*NA
        /nCarbonPAM_.value()*MCarbon)*pow(6.0*nCarbonPAM_.value()*MCarbon
        /(Foam::constant::mathematical::pi*rhoSoot_.value()
        *NA), 2.0/3.0)*NA*sqr(pamConcentration_[celli])
        *pow(abscissaNucleation, momentOrder);
}

// ************************************************************************* //
