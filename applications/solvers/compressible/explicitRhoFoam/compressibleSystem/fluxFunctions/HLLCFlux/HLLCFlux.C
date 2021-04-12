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

#include "HLLCFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxFunctions
{
    defineTypeNameAndDebug(HLLC, 0);
    addToRunTimeSelectionTable(fluxFunction, HLLC, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFunctions::HLLC::HLLC(const fvMesh& mesh)
:
    fluxFunction(mesh),
    thermo_
    (
        mesh.lookupObject<rhoThermo>("thermophysicalProperties")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxFunctions::HLLC::~HLLC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxFunctions::HLLC::updateFluxes
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

    volScalarField gamma("gamma", thermo_.gamma());

    surfaceScalarField gammaOwn
    (
        fvc::interpolate(gamma, own_, schemeName(gamma.name()))
    );

    surfaceScalarField gammaNei
    (
        fvc::interpolate(gamma, nei_, schemeName(gamma.name()))
    );

    surfaceScalarField aOwn
    (
        fvc::interpolate(a, own_, schemeName(a.name()))
    );

    surfaceScalarField aNei
    (
        fvc::interpolate(a, nei_, schemeName(a.name()))
    );

    surfaceScalarField UvOwn(UOwn & normal);
    UvOwn.setOriented(false);

    surfaceScalarField UvNei(UNei & normal);
    UvNei.setOriented(false);

    surfaceScalarField EOwn("EOwn", HOwn - pOwn/rhoOwn);
    surfaceScalarField ENei("ENei", HNei - pNei/rhoNei);

    // Averages
    surfaceScalarField aBar("aBar", 0.5*(aOwn + aNei));
    surfaceScalarField rhoBar("rhoBar", 0.5*(rhoOwn + rhoNei));

    // Estimate pressure
    surfaceScalarField pStar
    (
        "pStar",
        max
        (
            0.5*(pOwn + pNei - rhoBar*aBar*(UvNei - UvOwn)),
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );

    // Compute wave speeds
    surfaceScalarField rOwn
    (
        "rOwn",
        pos0(pOwn - pStar)
      + neg(pOwn - pStar)
       *sqrt(1.0 + (gammaOwn + 1.0)/(2.0*gammaOwn)*(pStar/pOwn - 1.0))
    );

    surfaceScalarField rNei
    (
        "rNei",
        pos0(pNei - pStar)
      + neg(pNei - pStar)
       *sqrt(1.0 + (gammaNei + 1.0)/(2.0*gammaNei)*(pStar/pNei - 1.0))
    );

    surfaceScalarField SOwn("SOwn", UvOwn - aOwn*rOwn);
    surfaceScalarField SNei("SNei", UvNei + aNei*rNei);

    surfaceScalarField SStar
    (
        "SStar",
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )/(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei))
    );

    // Compute contact pressure
    surfaceScalarField pOwnNei
    (
        "pOwnNei",
        0.5
       *(
            pOwn + pNei
          + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn)
          + rhoNei*(SNei - UvNei)*(SStar - UvNei)
        )
    );

    // Weighting factors for computing right and left star states
    surfaceScalarField rDeltaSOwn("rDeltaSOwn", 1/(SOwn - SStar));
    surfaceScalarField rDeltaSNei("rDeltaSNei", 1/(SNei - SStar));

    // Compute fluxes
    // Mass
    surfaceScalarField massFluxOwn(rhoOwn*UvOwn);

    surfaceScalarField massFluxStarOwn
    (
        SStar*(SOwn*rhoOwn - massFluxOwn)*rDeltaSOwn
    );

    surfaceScalarField massFluxNei(rhoNei*UvNei);

    surfaceScalarField massFluxStarNei
    (
        SStar*(SNei*rhoNei - massFluxNei)*rDeltaSNei
    );

    // Momentum
    surfaceVectorField pOwnTimesNormal(pOwn*normal);
    pOwnTimesNormal.setOriented(false);

    surfaceVectorField pOwnNeiTimesNormal(pOwnNei*normal);
    pOwnNeiTimesNormal.setOriented(false);

    surfaceVectorField momentumFluxOwn(UOwn*massFluxOwn + pOwnTimesNormal);

    surfaceVectorField momentumFluxStarOwn
    (
        (
            SStar*(SOwn*rhoOwn*UOwn - momentumFluxOwn)
          + SOwn*pOwnNeiTimesNormal
        )*rDeltaSOwn
    );

    surfaceVectorField pNeiTimesNormal(pNei*normal);
    pNeiTimesNormal.setOriented(false);

    surfaceVectorField momentumFluxNei(UNei*massFluxNei + pNeiTimesNormal);

    surfaceVectorField momentumFluxStarNei
    (
        (
            SStar*(SNei*rhoNei*UNei - momentumFluxNei)
          + SNei*pOwnNeiTimesNormal
        )*rDeltaSNei
    );

    // Energy
    surfaceScalarField energyFluxOwn(HOwn*massFluxOwn);

    surfaceScalarField energyFluxStarOwn
    (
        (
            SStar*(SOwn*rhoOwn*EOwn - energyFluxOwn)
          + SOwn*pOwnNei*SStar
        )*rDeltaSOwn
    );

    surfaceScalarField energyFluxNei(HNei*massFluxNei);
    
    surfaceScalarField energyFluxStarNei
    (
        (
            SStar*(SNei*rhoNei*ENei - energyFluxNei)
          + SNei*pOwnNei*SStar
        )*rDeltaSNei
    );

    // Set total fluxes
    massFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*massFluxOwn
          + pos0(SStar)*neg(SOwn)*massFluxStarOwn
          + pos0(SNei)*neg(SStar)*massFluxStarNei
          + neg(SNei)*massFluxNei
        );

    massFlux.setOriented(true);

    momentumFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*momentumFluxOwn
          + pos0(SStar)*neg(SOwn)*momentumFluxStarOwn
          + pos0(SNei)*neg(SStar)*momentumFluxStarNei
          + neg(SNei)*momentumFluxNei
        );

    momentumFlux.setOriented(true);

    energyFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*energyFluxOwn
          + pos0(SStar)*neg(SOwn)*energyFluxStarOwn
          + pos0(SNei)*neg(SStar)*energyFluxStarNei
          + neg(SNei)*energyFluxNei
        );

    energyFlux.setOriented(true);
}
