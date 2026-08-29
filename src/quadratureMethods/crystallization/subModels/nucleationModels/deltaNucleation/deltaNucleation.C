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
    Copyright (C) 2019-2023 Alberto Passalacqua
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

#include "deltaNucleation.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(deltaNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        deltaNucleation,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::nucleationModels::deltaNucleation::deltaNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    kb_("kb", inv(dimVolume*dimTime), dict),
    powCoeff_(readScalar(dict.lookup("powCoeff"))),
    sigma_
    (
        mesh.lookupObject<volScalarField>("sigma")
    ),
    JField_
    (
        IOobject
        (
            "nucleationRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", inv(dimVolume*dimTime), 0.0)
    ),
    d_nucleation_
    (
        "d_nucleation",
        dimLength,
        dict.lookupOrDefault<scalar>("d_nucleation", dict.lookupOrDefault<scalar>("deltaWidth", 1.0e-6))
    )

{}

Foam::populationBalanceSubModels::nucleationModels::deltaNucleation::~deltaNucleation()
{}

Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::deltaNucleation::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    
    const scalar sigmaVal = sigma_[celli];

    const scalar d_nucleation = d_nucleation_.value();

    scalar J = 0.0;
    if (sigmaVal > 0.0)
    {
        J = kb_.value() * std::pow(sigmaVal, powCoeff_);
    }
    
    const scalar nucleationSource = J * std::pow(d_nucleation, momentOrder);

    const_cast<volScalarField&>(JField_)[celli] = J;

    static label debugCounter = 0;
    if (debug && celli < 3 && debugCounter % 1000 == 0)
    {
        Pout<< "Cell " << celli << ": sigma=" << sigmaVal 
            << ", J=" << J << " #/(m³·s)"
            << ", moment_" << momentOrder << " source=" << nucleationSource << endl;
    }
    if (celli == 0) debugCounter++;

    return nucleationSource;
}

// ************************************************************************* //
