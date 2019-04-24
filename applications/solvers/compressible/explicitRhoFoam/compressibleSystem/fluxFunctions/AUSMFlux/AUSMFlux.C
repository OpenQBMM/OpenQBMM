/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Jeff Heylmun
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

#include "AUSMFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxFunctions
{
    defineTypeNameAndDebug(AUSM, 0);
    addToRunTimeSelectionTable(fluxFunction, AUSM, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFunctions::AUSM::AUSM(const fvMesh& mesh)
:
    fluxFunction(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxFunctions::AUSM::~AUSM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxFunctions::AUSM::updateFluxes
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
        fvc::interpolate(rho, own_, interpScheme(rho.name()))
    );
    surfaceScalarField rhoNei
    (
        fvc::interpolate(rho, nei_, interpScheme(rho.name()))
    );

    surfaceVectorField UOwn(fvc::interpolate(U, own_, interpScheme(U.name())));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, interpScheme(U.name())));

    surfaceScalarField HOwn(fvc::interpolate(H, own_, interpScheme(H.name())));
    surfaceScalarField HNei(fvc::interpolate(H, nei_, interpScheme(H.name())));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, interpScheme(p.name())));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, interpScheme(p.name())));

    surfaceScalarField aOwn
    (
        fvc::interpolate(a, own_, interpScheme(a.name()))
    );
    surfaceScalarField aNei
    (
        fvc::interpolate(a, nei_, interpScheme(a.name()))
    );

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    // Compute slpit Mach numbers
    surfaceScalarField MaOwn("MaOwn", UvOwn/max(aOwn, minU));
    surfaceScalarField MaNei("MaNei", UvNei/max(aNei, minU));
    surfaceScalarField magMaOwn(mag(MaOwn));
    surfaceScalarField magMaNei(mag(MaNei));

    surfaceScalarField MaPlus
    (
        "MaPlus",
        neg(magMaOwn - 1.0)*(0.25*sqr(MaOwn + 1.0) + 0.125*(sqr(MaOwn) - 1.0))
      + pos0(MaOwn - 1.0)*MaOwn
    );
    surfaceScalarField MaMinus
    (
        "MaMinus",
        pos0(-1.0 - MaNei)*MaNei
      - neg(magMaNei - 1.0)*(0.25*sqr(MaNei - 1.0) - 0.125*(sqr(MaNei) - 1.0))
    );
    surfaceScalarField Ma12
    (
        "Ma12",
        MaPlus + MaMinus
    );

    // Compute slpit pressures
    surfaceScalarField pPlus
    (
        "pPlus",
        neg(magMaOwn - 1.0)*0.25*sqr(MaOwn + 1.0)*(2.0 - MaOwn)
      + pos0(MaOwn - 1.0)
    );
    surfaceScalarField pMinus
    (
        "pMinus",
        pos0(-1.0 - MaNei)
      + neg(magMaNei - 1.0)*0.25*sqr(MaNei - 1.0)*(2.0 + MaNei)
    );
    surfaceScalarField p12
    (
        "p12",
        pOwn*pPlus + pNei*pMinus
    );

    surfaceScalarField plus("plus", pos0(Ma12));
    surfaceScalarField minus("minus", neg(Ma12));

    massFlux =
        mesh_.magSf()
       *(
           max(0.0, Ma12)*rhoOwn*aOwn
         + min(0.0, Ma12)*rhoNei*aNei
        );

    momentumFlux =
        mesh_.magSf()
       *(
            max(0.0, Ma12)*rhoOwn*aOwn*UOwn
          + min(0.0, Ma12)*rhoNei*aNei*UNei
        )
      + p12*mesh_.Sf();

    energyFlux =
        mesh_.magSf()
       *(
           max(0.0, Ma12)*rhoOwn*aOwn*HOwn
         + min(0.0, Ma12)*rhoNei*aNei*HNei
        );
}
