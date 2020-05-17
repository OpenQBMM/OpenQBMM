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

#include "mixingKernel.H"

#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
    defineTypeNameAndDebug(mixingKernel, 0);
    defineRunTimeSelectionTable(mixingKernel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernel::mixingKernel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    Cphi_
    (
        dict.lookupOrDefault
        (
            "Cphi",
            dimensionedScalar("CPhiDefault", dimless, 2.0)
        )
    ),
    Cmixing_
    (
        dict.lookupOrDefault
        (
            "Cmixing",
            dimensionedScalar("CmixingDefault", dimless, 1.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernel::~mixingKernel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::mixingSubModels::mixingKernel::k() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if 
    (
        mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName)
    )
    {
        const cmpTurbModel& turb =
            mesh_.lookupObject<cmpTurbModel>
            (
                cmpTurbModel::propertiesName
            );

        return turb.k();
    }
    else if
    (
        mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName)
    )
    {
        const incompressible::turbulenceModel& turb =
            mesh_.lookupObject<icoTurbModel>
            (
                icoTurbModel::propertiesName
            );

        return turb.k();
    }

    FatalErrorInFunction
            << "No valid turbulence model found."
            << exit(FatalError);

    return volScalarField::null();
}

Foam::tmp<Foam::volScalarField> 
Foam::mixingSubModels::mixingKernel::epsilon() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if 
    (
        mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName)
    )
    {
        const cmpTurbModel& turb =
            mesh_.lookupObject<cmpTurbModel>
            (
                cmpTurbModel::propertiesName
            );

        return turb.epsilon();
    }
    else if
    (
        mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName)
    )
    {
        const incompressible::turbulenceModel& turb =
            mesh_.lookupObject<icoTurbModel>
            (
                icoTurbModel::propertiesName
            );

        return turb.epsilon();
    }

    FatalErrorInFunction
            << "No valid turbulence model found."
            << exit(FatalError);

    return volScalarField::null();
}

// ************************************************************************* //
