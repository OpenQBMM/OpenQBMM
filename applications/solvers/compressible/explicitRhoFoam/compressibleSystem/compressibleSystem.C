/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Alberto Passalacqua
     \\/     M anipulation  |
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

#include "compressibleSystem.H"
#include "fluxIntegrator.H"
#include "fluxFunction.H"
#include "surfaceInterpolate.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "constants.H"
#include "fvc.H"
#include "fvm.H"


// * * * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleSystem::compressibleSystem
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    thermoPtr_(rhoThermo::New(mesh)),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermoPtr_().rho()
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_(thermoPtr_->p()),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermoPtr_->he() + 0.5*magSqr(U_)
    ),
    H_
    (
        IOobject
        (
            "H",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        E_ + p_/rho_
    ),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U_
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*E_
    ),
    massFlux_
    (
        IOobject
        (
            "massFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimVelocity*dimDensity*dimArea, 0)
    ),
    momentumFlux_
    (
        IOobject
        (
            "momentumFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0",
            sqr(dimVelocity)*dimDensity*dimArea,
            Zero
        )
    ),
    energyFlux_
    (
        IOobject
        (
            "energyFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", pow3(dimVelocity)*dimDensity*dimArea, 0)
    )
{
    const word phiName = "phi";

    IOobject phiHeader
    (
        phiName,
        mesh_.time().timeName(),
        mesh_,
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
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
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
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U_),
                phiTypes
            )
        );
    }
    thermoPtr_->validate("compressibleSystem ", "e");
    encode();

    integrator_.set(new fluxIntegrator(*this));
    fluxFunction_ = fluxFunction::New(mesh_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleSystem::~compressibleSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibleSystem::setNSteps
(
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    label nSteps = storeFields.size();
    storedFieldIndexes_ = labelList(nSteps, label(-1));
    storedDeltaIndexes_ = labelList(nSteps, label(-1));

    label currFieldIndex = 0;
    label currDeltaIndex = 0;
    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        if (storeFields[stepi])
        {
            storedFieldIndexes_[stepi] = currFieldIndex;
            currFieldIndex++;

            rhos_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "rho" + Foam::name(stepi),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("zero", rho_.dimensions(), 0.0)
                )
            );
            rhoUs_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        "rhoU" + Foam::name(stepi),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedVector("zero", rhoU_.dimensions(), Zero)
                )
            );
            rhoEs_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "rhoE" + Foam::name(stepi),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("zero", rhoE_.dimensions(), 0.0)
                )
            );

        }

        if (storeDeltas[stepi])
        {
            storedDeltaIndexes_[stepi] = currDeltaIndex;
            currDeltaIndex++;

            deltaRhos_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "deltaRho" + Foam::name(stepi),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("zero", rho_.dimensions()/dimTime, 0.0)
                )
            );
            deltaRhoUs_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        "deltaRhoU" + Foam::name(stepi),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedVector("zero", rhoU_.dimensions()/dimTime, Zero)
                )
            );
            deltaRhoEs_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "deltaRhoE" + Foam::name(stepi),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("zero", rhoE_.dimensions()/dimTime, 0.0)
                )
            );

        }
    }
}

Foam::tmp<Foam::volScalarField>
Foam::compressibleSystem::speedOfSound() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "c",
            sqrt(thermoPtr_->gamma()*p_/rho_)
        )
    );
}

void Foam::compressibleSystem::advect
(
    const label stepi,
    const scalarList& coeffs,
    const scalarList& Fcoeffs,
    const dimensionedScalar& deltaT,
    const dimensionedVector& g
)
{
    label fi = storedFieldIndexes_[stepi];
    label di = storedDeltaIndexes_[stepi];

    // Store current step if needed later
    if (fi != -1)
    {
        rhos_[fi] = rho_;
        rhoUs_[fi] = rhoU_;
        rhoEs_[fi] = rhoE_;
    }

    volScalarField deltaRho(-fvc::div(massFlux_));
    volVectorField deltaRhoU(-fvc::div(momentumFlux_) + rho_*g);
    volScalarField deltaRhoE(-fvc::div(energyFlux_) + (rho_*g & U_));

    // Store deltas if needed later
    if (di != -1)
    {
        deltaRhos_[di] = deltaRho;
        deltaRhoUs_[di] = deltaRhoU;
        deltaRhoEs_[di] = deltaRhoE;
    }

    volScalarField rho(rho_*coeffs[stepi]);
    volVectorField rhoU(rhoU_*coeffs[stepi]);
    volScalarField rhoE(rhoE_*coeffs[stepi]);

    deltaRho *= Fcoeffs[stepi];
    deltaRhoU *= Fcoeffs[stepi];
    deltaRhoE *= Fcoeffs[stepi];

    label fieldi = 0;
    label deltai = 0;
    for (label i = 0; i < stepi; i++)
    {
        if (storedFieldIndexes_[i] != -1)
        {
            rho += rhos_[fieldi]*coeffs[i];
            rhoU += rhoUs_[fieldi]*coeffs[i];
            rhoE += rhoEs_[fieldi]*coeffs[i];
            fieldi++;
        }

        if (storedDeltaIndexes_[i] != -1)
        {
            deltaRho += deltaRhos_[deltai]*Fcoeffs[i];
            deltaRhoU += deltaRhoUs_[deltai]*Fcoeffs[i];
            deltaRhoE += deltaRhoEs_[deltai]*Fcoeffs[i];
            deltai++;
        }
    }

    rho_ = rho + deltaT*deltaRho;
    rho_.correctBoundaryConditions();

    rhoU_ = rhoU + deltaT*deltaRhoU;
    //- Ensure unused directions are zero
    rhoU_ =
        cmptMultiply
        (
            rhoU_,
            (vector(mesh_.solutionD()) + vector::one)/2.0
        );
    rhoU_.correctBoundaryConditions();

    rhoE_ = rhoE + deltaT*deltaRhoE;
    rhoE_.correctBoundaryConditions();
}

void Foam::compressibleSystem::integrateFluxes
(
    const dimensionedVector& g
)
{
    integrator_->integrateFluxes(g);
}


void Foam::compressibleSystem::updateFluxes()
{

    // calculate fluxes with
    fluxFunction_->updateFluxes
    (
        massFlux_,
        momentumFlux_,
        energyFlux_,
        rho_,
        U_,
        H_,
        p_,
        speedOfSound()()
    );
}


void Foam::compressibleSystem::decode()
{
    thermoPtr_->rho() = rho_;

    U_ = rhoU_/rho_;
    U_.correctBoundaryConditions();

    phi() = fvc::flux(U_);

    E_ = rhoE_/rho_;
    thermoPtr_->he() = E_ - 0.5*magSqr(U_);

    thermoPtr_->correct();
    p_ = rho_/thermoPtr_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermoPtr_->psi().boundaryField()*p_.boundaryField();

    H_ = E_ + p_/rho_;
}


void Foam::compressibleSystem::encode()
{
    rho_ = thermoPtr_->rho();
    rho_.correctBoundaryConditions();

    rhoU_ = rho_*U_;
    rhoU_.boundaryFieldRef() == rho_.boundaryField()*U_.boundaryField();

    rhoE_ = rho_*E_;
    rhoE_.boundaryFieldRef() ==
        rho_.boundaryField()*
        (
            E_.boundaryField() + 0.5*magSqr(U_.boundaryField())
        );
}


void Foam::compressibleSystem::correctThermo()
{
    E_ = thermoPtr_->he() + 0.5*magSqr(U_);

    thermoPtr_->correct();
    p_ = rho_/thermoPtr_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermoPtr_->psi().boundaryField()*p_.boundaryField();

    H_ = E_ + p_/rho_;
}


// ************************************************************************* //
