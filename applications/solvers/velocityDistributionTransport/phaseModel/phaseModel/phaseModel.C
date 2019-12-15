/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
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

#include "phaseModel.H"
#include "fvMatrix.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha", dimless, Zero)
    ),
    name_(phaseName),
    phaseDict_
    (
        dict.subDict(name_)
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        phaseDict_
    ),
    alphaMax_(phaseDict_.lookupOrDefault<scalar>("alphaMax", 1)),
    thermo_(rhoThermo::New(mesh, name_)),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity, Zero)
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("d", dimLength, phaseDict_)
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alphaPhiZero", dimensionSet(0, 3, -1, 0, 0), 0.0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alphaRhoPhiZero", dimensionSet(1, 0, -1, 0, 0), 0.0)
    )
{
    thermo_->validate("phaseModel " + name_, "h", "e");

    U_.correctBoundaryConditions();
    const word phiName = IOobject::groupName("phi", name_);
    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvsPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U_),
                phiTypes
            )
        );
    }

    turbulence_ =
        PhaseCompressibleTurbulenceModel<phaseModel>::New
        (
            *this,
            rho(),
            U_,
            alphaRhoPhi_,
            phiPtr_(),
            *this
        );
}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence()
{
    return turbulence_();
}


const Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence() const
{
    return turbulence_();
}


Foam::tmp<Foam::volVectorField> Foam::phaseModel::Vs(const label nodei) const
{
    return
        tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "zero",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh(),
                dimensionedVector("0", dimVelocity, Zero)
            )
        );
}

// ************************************************************************* //
