/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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

#include "firstOrderFluidVelocityAdvection.H"
#include "upwind.H"
#include "fvcSurfaceIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{
    defineTypeNameAndDebug(firstOrderFluidVelocity, 0);

    addToRunTimeSelectionTable
    (
        velocityMomentAdvection,
        firstOrderFluidVelocity,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityAdvection::firstOrderFluidVelocity::firstOrderFluidVelocity
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const List<supportType>& supports
)
:
    velocityMomentAdvection(dict, quadrature, supports),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    mesh_(quadrature.moments()(0).mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityAdvection::firstOrderFluidVelocity::~firstOrderFluidVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::surfaceScalarField&
Foam::velocityAdvection::firstOrderFluidVelocity::phi() const
{
    if (!mesh_.foundObject<surfaceScalarField>(phiName_))
    {
        FatalErrorInFunction
            << "The flux " << phiName_ << " the moments of " << name_
            << " are carried by is not in the registry of the mesh." << nl
            << "    Name the flux of the fluid with phi in the "
            << "momentAdvection dictionary of the quadrature." << nl
            << exit(FatalError);
    }

    return mesh_.lookupObject<surfaceScalarField>(phiName_);
}


bool Foam::velocityAdvection::firstOrderFluidVelocity::densityWeighted() const
{
    return phi().dimensions() == dimMass/dimTime;
}


Foam::scalar
Foam::velocityAdvection::firstOrderFluidVelocity::realizableCo() const
{
    return 1.0;
}


Foam::scalar Foam::velocityAdvection::firstOrderFluidVelocity::CoNum() const
{
    const surfaceScalarField& phiFluid = phi();

    const scalarField sumPhi
    (
        fvc::surfaceSum(mag(phiFluid))().primitiveField()
    );

    // The volume of fluid through the faces of a cell over its volume, or
    // the mass through them over its mass
    if (densityWeighted())
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        return
            0.5*gMax(sumPhi/(rho.primitiveField()*mesh_.V().field()))
           *mesh_.time().deltaTValue();
    }

    return
        0.5*gMax(sumPhi/mesh_.V().field())
       *mesh_.time().deltaTValue();
}


void Foam::velocityAdvection::firstOrderFluidVelocity::update()
{
    update(phi(), false);
}


void Foam::velocityAdvection::firstOrderFluidVelocity::update
(
    const surfaceScalarField& phiFluid,
    const bool wallCollisions
)
{
    // Every node moves with the fluid, so the moments of a face are those
    // of the cell upwind of it, and the advection term of every moment is
    // the divergence of its upwind flux
    const upwind<scalar> faceValues(mesh_, phiFluid);

    forAll(divMoments_, momenti)
    {
        // The divergence takes the dimensions of the flux: of a volume, or
        // of a mass
        divMoments_[momenti].dimensions().reset
        (
            moments_[momenti].dimensions()*phiFluid.dimensions()/dimVolume
        );

        divMoments_[momenti] =
            fvc::surfaceIntegrate(faceValues.flux(moments_[momenti]));
    }
}


void Foam::velocityAdvection::firstOrderFluidVelocity::update
(
    const mappedPtrList<volVectorField>& Us,
    const bool wallCollisions
)
{
    FatalErrorInFunction
        << "The nodes of this distribution share the velocity of the fluid, "
        << "so there is no velocity of a node to advect them with." << nl
        << exit(FatalError);
}

// ************************************************************************* //
