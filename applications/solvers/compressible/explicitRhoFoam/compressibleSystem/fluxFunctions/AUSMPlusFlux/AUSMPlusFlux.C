/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 Alberto Passalacqua
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

#include "AUSMPlusFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxFunctions
{
    defineTypeNameAndDebug(AUSMPlus, 0);
    addToRunTimeSelectionTable(fluxFunction, AUSMPlus, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFunctions::AUSMPlus::AUSMPlus(const fvMesh& mesh)
:
    fluxFunction(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxFunctions::AUSMPlus::~AUSMPlus()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxFunctions::AUSMPlus::updateFluxes
(
    surfaceScalarField& massFlux,
    surfaceVectorField& momentumFlux,
    surfaceScalarField& energyFlux,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& H,
    const volScalarField& p,
    const volScalarField& a
)
{
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());
    
    dimensionedScalar minU("smallU", dimVelocity, SMALL);

    surfaceScalarField rhoOwn
    (
        fvc::interpolate(rho, own_, schemeName(rho.name()))
    );

    surfaceScalarField rhoNei
    (
        fvc::interpolate(rho, nei_, schemeName(rho.name()))
    );

    surfaceVectorField UOwn(fvc::interpolate(U, own_, schemeName(U.name())));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, schemeName(U.name())));

    surfaceScalarField HOwn(fvc::interpolate(H, own_, schemeName(H.name())));
    surfaceScalarField HNei(fvc::interpolate(H, nei_, schemeName(H.name())));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, schemeName(p.name())));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, schemeName(p.name())));

    surfaceScalarField aOwn
    (
        fvc::interpolate(a, own_, schemeName(a.name()))
    );

    surfaceScalarField aNei
    (
        fvc::interpolate(a, nei_, schemeName(a.name()))
    );

    surfaceScalarField UvOwn(UOwn & normal);
    UvOwn.setOriented(false); // Specific to OF+

    surfaceScalarField UvNei(UNei & normal);
    UvNei.setOriented(false); // Specific to OF+

    // Compute split Mach numbers
    surfaceScalarField MaOwn("MaOwn", UvOwn/max(aOwn, minU));
    surfaceScalarField MaNei("MaNei", UvNei/max(aNei, minU));
    surfaceScalarField magMaOwn(mag(MaOwn));
    surfaceScalarField magMaNei(mag(MaNei));

    surfaceScalarField deltapOwn
    (
        "deltapOwn",
        pos0(magMaOwn - 1)*sign(MaOwn)
      + neg(magMaOwn - 1)
       *MaOwn/2.0
       *(
            3.0
          - sqr(MaOwn)
          + 4.0*alpha_*sqr(sqr(MaOwn - 1.0))
        )
    );

    surfaceScalarField deltapNei
    (
        "deltapNei",
        pos0(magMaNei - 1)*sign(MaNei)
      + neg(magMaNei - 1)
       *MaNei/2.0
       *(
            3.0
          - sqr(MaNei)
          + 4.0*alpha_*sqr(sqr(MaNei - 1.0))
        )
    );

    surfaceScalarField deltap("deltap", deltapOwn*pOwn + deltapNei*pNei);

    surfaceScalarField magMachOwn
    (
        "magMachOwn",
        pos0(magMaOwn - 1)*magMaOwn
      + neg(magMaOwn - 1)
       *(
            0.5*(sqr(MaOwn) + 1.0)
          + 2.0*beta_*sqr(sqr(MaOwn) - 1.0)
        )
    );

    surfaceScalarField magMachNei
    (
        "magMachNei",
        pos0(magMaNei - 1)*magMaNei
      + neg(magMaNei - 1)
       *(
            0.5*(sqr(MaNei) + 1.0)
          + 2.0*beta_*sqr(sqr(MaNei) - 1.0)
        )
    );

    surfaceScalarField deltaMa12
    (
        "deltaMa12",
        magMachOwn - magMachNei
    );

    surfaceScalarField Ma12
    (
        "Ma12",
        MaOwn + MaNei - deltaMa12
    );

    surfaceScalarField a12("a12", sqrt(aOwn*aNei));
    surfaceScalarField rhoPhi(fvc::interpolate(rho*U) & normal);

    surfaceVectorField rhoUPhi
    (
        (fvc::interpolate(rho*U*U) & normal)
      + fvc::interpolate(p)*normal
    );
    
    surfaceScalarField rhoHPhi(fvc::interpolate(rho*H*U) & normal);
    
    surfaceScalarField a12DeltaMaRho
    (
       0.5*a12
      *(
           (0.5*deltaMa12 - mag(Ma12))*rhoOwn
         + (0.5*deltaMa12 + mag(Ma12))*rhoNei
       )
    );

    a12DeltaMaRho.setOriented(true);

    massFlux =
        mesh_.magSf()
       *(
            rhoPhi
          - a12DeltaMaRho
        );

    surfaceVectorField a12DeltaMaRhoU
    (
        0.5*a12
       *(
            (0.5*deltaMa12 - mag(Ma12))*rhoOwn*UOwn
          + (0.5*deltaMa12 + mag(Ma12))*rhoNei*UNei
        )
    );

    a12DeltaMaRhoU.setOriented(true);

    momentumFlux =
        mesh_.magSf()
       *(
            rhoUPhi
          - a12DeltaMaRhoU
        );
      - 0.5*deltap*mesh_.Sf();

    surfaceScalarField a12DeltaMaRhoH
    (
        0.5*a12
      *(
           (0.5*deltaMa12 - mag(Ma12))*rhoOwn*HOwn
         + (0.5*deltaMa12 + mag(Ma12))*rhoNei*HNei
       )
    );

    a12DeltaMaRhoH.setOriented(true);
    
    energyFlux =
        mesh_.magSf()
       *(
            rhoHPhi
          - a12DeltaMaRhoH
        );
}
