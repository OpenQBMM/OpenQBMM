/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added additional return functions so that class can
                            be extended to polydisperse
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
#include "twoPhaseSystem.H"
#include "fvMatrix.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "dragModel.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& dict,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("alpha", dimless, 0)
    ),
    fluid_(fluid),
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
    thermo_(rhoThermo::New(fluid.mesh(), name_)),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName
            (
                "d",
                name_
            ),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar
        (
            "d",
            dimLength,
            phaseDict_.lookupOrDefault<scalar>
            (
                "d",
                0.0
            )
        )
    ),
    BGviscosity_(phaseDict_.lookupOrDefault("BGviscosity", false))
{
    thermo_->validate("phaseModel " + name_, "h", "e");

    const word phiName = IOobject::groupName("phi", name_);

    IOobject phiHeader
    (
        phiName,
        fluid.mesh().time().timeName(),
        fluid.mesh(),
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
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh()
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
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
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
            thermo_->rho(),
            U_,
            alphaRhoPhi_,
            phi(),
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

const Foam::phaseModel& Foam::phaseModel::otherPhase() const
{
    return fluid_.otherPhase(*this);
}


void Foam::phaseModel::correct()
{
    return;
}


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
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fluid_.mesh(),
                dimensionedVector("0", dimVelocity, Zero)
            )
        );
}


void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();

    // Ensure that the flux at inflow BCs is preserved
    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            const scalarField& phip = phi().boundaryField()[patchi];
            const scalarField& alphap = boundaryField()[patchi];

            forAll(alphaPhip, facei)
            {
                if (phip[facei] < SMALL)
                {
                    alphaPhip[facei] = alphap[facei]*phip[facei];
                }
            }
        }
    }
}


// ************************************************************************* //
