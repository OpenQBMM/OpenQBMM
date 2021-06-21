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

#include "HLLFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxFunctions
{
    defineTypeNameAndDebug(HLL, 0);
    addToRunTimeSelectionTable(fluxFunction, HLL, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFunctions::HLL::HLL(const fvMesh& mesh)
:
    fluxFunction(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxFunctions::HLL::~HLL()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxFunctions::HLL::updateFluxes
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

    surfaceScalarField aOwn(fvc::interpolate(a, own_, schemeName(a.name())));
    surfaceScalarField aNei(fvc::interpolate(a, nei_, schemeName(a.name())));

    surfaceScalarField UvOwn(UOwn & normal);
    UvOwn.setOriented(false);

    surfaceScalarField UvNei(UNei & normal);
    UvNei.setOriented(false);

    surfaceScalarField EOwn("EOwn", HOwn - pOwn/rhoOwn);
    surfaceScalarField ENei("ENei", HNei - pNei/rhoNei);

    // Averages
    surfaceScalarField aBar("aBar", 0.5*(aOwn + aNei));
    surfaceScalarField rhoBar("rhoBar", 0.5*(rhoOwn + rhoNei));

    // Fluxes
    surfaceScalarField massFluxOwn(rhoOwn*UvOwn);
    surfaceScalarField massFluxNei(rhoNei*UvNei);

    surfaceVectorField pOwnTimesNormal(pOwn*normal);
    pOwnTimesNormal.setOriented(false);

    surfaceVectorField pNeiTimesNormal(pNei*normal);
    pNeiTimesNormal.setOriented(false);

    surfaceVectorField momentumFluxOwn(UOwn*massFluxOwn + pOwnTimesNormal);
    surfaceVectorField momentumFluxNei(UNei*massFluxNei + pNeiTimesNormal);

    surfaceScalarField energyFluxOwn(HOwn*massFluxOwn);
    surfaceScalarField energyFluxNei(HNei*massFluxNei);

    surfaceScalarField SOwn("SOwn", min(UvOwn - aOwn, UvNei - aNei));
    surfaceScalarField SNei("SNei", max(UvNei + aNei, UvOwn + aNei));

    surfaceScalarField rDeltaS("rDeltaS", 1.0/(SNei - SOwn));

    // Compute fluxes
    massFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*massFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*massFluxOwn - SOwn*massFluxNei
              + SOwn*SNei*(rhoNei - rhoOwn)
            )*rDeltaS
          + neg(SNei)*massFluxNei
        );

    momentumFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*momentumFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*momentumFluxOwn - SOwn*momentumFluxNei
              + SOwn*SNei*(rhoNei*UNei - rhoOwn*UOwn)
            )*rDeltaS
          + neg(SNei)*momentumFluxNei
        );

    energyFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*energyFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*energyFluxOwn - SOwn*energyFluxNei
              + SOwn*SNei*(rhoNei*ENei - rhoOwn*EOwn)
            )*rDeltaS
          + neg(SNei)*energyFluxNei
        );
}
