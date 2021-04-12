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

#include "RoeFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxFunctions
{
    defineTypeNameAndDebug(Roe, 0);
    addToRunTimeSelectionTable(fluxFunction, Roe, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFunctions::Roe::Roe(const fvMesh& mesh) : fluxFunction(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxFunctions::Roe::~Roe()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxFunctions::Roe::updateFluxes
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

    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());

    surfaceScalarField wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    surfaceScalarField wNei(1.0 - wOwn);

    surfaceScalarField rhoTilde
    (
        "rhoTilde",
        sqrt
        (
            max
            (
                rhoOwn*rhoNei,
                dimensionedScalar("SMALL", sqr(dimDensity), SMALL)
            )
        )
    );

    surfaceVectorField UTilde
    (
        "UTilde",
        UOwn*wOwn + UNei*wNei
    );

    UTilde.setOriented(true);

    surfaceScalarField UvTilde(UTilde & normal);
    UvTilde.setOriented(false);

    surfaceScalarField HTilde
    (
        "HTilde",
        HOwn*wOwn + HNei*wNei
    );

    surfaceScalarField aTilde
    (
        "aTilde",
        aOwn*wOwn + aNei*wNei
    );

    surfaceScalarField deltaRho(rhoNei - rhoOwn);
    surfaceVectorField deltaU(UNei - UOwn);

    surfaceScalarField deltaUv(deltaU & normal);
    deltaUv.setOriented(false);

    surfaceScalarField deltaP(pNei - pOwn);

    surfaceScalarField lambda1(mag(UvTilde));
    surfaceScalarField lambda2(mag(UvTilde + aTilde));
    surfaceScalarField lambda3(mag(UvTilde - aTilde));

    surfaceScalarField alpha1
    (
        deltaRho - deltaP/sqr(aTilde)
    );

    surfaceScalarField alpha2
    (
        (deltaP + rhoTilde*aTilde*deltaUv)/(2.0*sqr(aTilde))
    );
    
    surfaceScalarField alpha3
    (
        (deltaP - rhoTilde*aTilde*deltaUv)/(2.0*sqr(aTilde))
    );

    // U Row
    surfaceVectorField K21(UTilde);
    surfaceVectorField K224((UTilde + aTilde*normal));
    surfaceVectorField K25((UTilde - aTilde*normal));

    // E row
    surfaceScalarField K31(0.5*magSqr(UTilde));
    surfaceScalarField K324((HTilde + aTilde*UvTilde));
    surfaceScalarField K35((HTilde - aTilde*UvTilde));

    // Compute owner and neighbour fluxes
    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

//     forAll(lambda1, facei)
//     {
//         scalar eps = 0.1*aTilde[facei]; //adjustable parameter
//
//         if
//         (
//             lambda1[facei] < eps
//          || lambda2[facei] < eps
//          || lambda3[facei] < eps
//         )
//         {
//             lambda1[facei] = (sqr(lambda1[facei]) + sqr(eps))/(2.0*eps);
//             lambda2[facei] = (sqr(lambda2[facei]) + sqr(eps))/(2.0*eps);
//             lambda3[facei] = (sqr(lambda3[facei]) + sqr(eps))/(2.0*eps);
//         }
//
//         // First eigenvalue: U - c
//         eps =
//             2.0*max(0,(UvNei[facei] - aNei[facei])
//           - (UvOwn[facei] - aOwn[facei]));
//         if (lambda1[facei] < eps)
//         {
//             lambda1[facei] = (sqr(lambda1[facei]) + sqr(eps))/(2.0*eps);
//         }
//
//         // Second eigenvalue: U
//         eps = 2*max(0, UvNei[facei] - UvOwn[facei]);
//         if (lambda2[facei] < eps)
//         {
//             lambda2[facei] = (sqr(lambda2[facei]) + sqr(eps))/(2.0*eps);
//         }
//
//         // Third eigenvalue: U + c
//         eps =
//             2.0*max(0,(UvNei[facei] + aNei[facei])
//           - (UvOwn[facei] + aOwn[facei]));
//         if (lambda3[facei] < eps)
//         {
//             lambda3[facei] = (sqr(lambda3[facei]) + sqr(eps))/(2.0*eps);
//         }
//     }

    surfaceScalarField massFluxOwn(rhoOwn*UvOwn);
    surfaceScalarField massFluxNei(rhoNei*UvNei);

    surfaceVectorField momentumFluxOwn
    (
        UOwn*massFluxOwn + pOwn*normal
    );

    surfaceVectorField momentumFluxNei
    (
        UNei*massFluxNei + pNei*normal
    );

    surfaceScalarField energyFluxOwn(HOwn*massFluxOwn);
    surfaceScalarField energyFluxNei(HNei*massFluxNei);

    surfaceScalarField lambdaAlpha
    (
        mag(lambda1)*alpha1
      + mag(lambda2)*alpha2
      + mag(lambda3)*alpha3
    );

    lambdaAlpha.setOriented(true);

    // Compute fluxes
    massFlux = 0.5*mesh_.magSf()*(massFluxOwn + massFluxNei - lambdaAlpha);
    massFlux.setOriented(true);

    surfaceVectorField lambdaAlphaK
    (
        mag(lambda1)*alpha1*K21
      + mag(lambda2)*alpha2*K224
      + mag(lambda3)*alpha3*K25
    );

    lambdaAlphaK.setOriented(true);

    momentumFlux =
        0.5*mesh_.magSf()*(momentumFluxOwn + momentumFluxNei - lambdaAlphaK);

    momentumFlux.setOriented(true);

    surfaceScalarField lambdaAlphaKE
    (
        mag(lambda1)*alpha1*K31
      + mag(lambda2)*alpha2*K324
      + mag(lambda3)*alpha3*K35
    );

    lambdaAlphaKE.setOriented(true);

    energyFlux =
        0.5*mesh_.magSf()*(energyFluxOwn + energyFluxNei - lambdaAlphaKE);

    energyFlux.setOriented(true);
}
