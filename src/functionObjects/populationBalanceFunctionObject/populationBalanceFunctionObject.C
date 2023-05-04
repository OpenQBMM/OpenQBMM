/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2023 Alberto Passalacqua
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

#include "populationBalanceFunctionObject.H"
#include "surfaceFields.H"
#include "fvcFlux.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(populationBalance, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        populationBalance,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::populationBalance::updateVolumetricFaceFlux() 
{
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVolume/dimTime)
    {
        phi_ = phi;
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        phi_ = phi/fvc::interpolate(rho);
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::populationBalance::populationBalance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    populationBalanceProperties_
    (
        IOobject
        (
            "populationBalanceProperties",
            runTime.constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    phi_
    (
        IOobject
        (
            "phiPopulationBalance",
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zeroPhi", dimVolume/dimTime, 0.0)
    )
{
    read(dict);

    updateVolumetricFaceFlux();

    populationBalanceModelPtr_ = populationBalanceModel::New
    (
        "populationBalance",
        populationBalanceProperties_,
        phi_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::populationBalance::~populationBalance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::populationBalance::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    return true;
}


Foam::wordList Foam::functionObjects::populationBalance::fields() const
{
    return wordList{phiName_};
}


bool Foam::functionObjects::populationBalance::execute()
{
    Info<< type() << " write:" << endl;

    // Update flux (face velocity)
    updateVolumetricFaceFlux();

    // Solve population balance
    populationBalanceModelPtr_->solve();

    Info<< endl;

    return true;
}

bool Foam::functionObjects::populationBalance::write()
{
    return true;
}


// ************************************************************************* //