/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-06-05  Jeff Heylmun:   Modified to allow for use of anisotropic Gaussian
                            model.
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

#include "equilibrium.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(equilibrium, 0);
    addToRunTimeSelectionTable
    (
        kineticTheoryModel,
        equilibrium,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::equilibrium::equilibrium
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    kineticTheoryModel(dict, phase)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::equilibrium::~equilibrium()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModels::equilibrium::solve
(
    const volScalarField& beta,
    const volScalarField& alpha,
    const volTensorField& gradU,
    const volSymmTensorField D
)
{
    // Local references
    const volScalarField& rho = phase_.rho();
    const volScalarField& da = phase_.d();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    // Equilibrium => dissipation == production
    // Eq. 4.14, p.82
    volScalarField K1("K1", 2.0*(1.0 + e_)*rho*g0_);
    volScalarField K3
    (
        "K3",
        0.5*da*rho*
        (
            (sqrtPi/(3.0*(3.0 - e_)))
            *(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha*g0_)
            +1.6*alpha*g0_*(1.0 + e_)/sqrtPi
        )
    );

    volScalarField K2
    (
        "K2",
        4.0*da*rho*(1.0 + e_)*alpha*g0_/(3.0*sqrtPi) - 2.0*K3/3.0
    );

    volScalarField K4("K4", 12.0*(1.0 - sqr(e_))*rho*g0_/(da*sqrtPi));

    volScalarField trD
    (
        "trD",
        alpha/(alpha + residualAlpha_)
       *fvc::div(phase_.phi()*this->h2f())
    );
    volScalarField tr2D("tr2D", sqr(trD));
    volScalarField trD2("trD2", tr(D & D));

    volScalarField t1("t1", K1*alpha + rho);
    volScalarField l1("l1", -t1*trD);
    volScalarField l2("l2", sqr(t1)*tr2D);
    volScalarField l3
    (
        "l3",
        4.0
        *K4
        *alpha
        *(2.0*K3*trD2 + K2*tr2D)
    );

    Theta_ +=
        sqr
        (
            (l1 + sqrt(l2 + l3))
            /(2.0*max(alpha, residualAlpha_)*K4)
        );


    Theta_.max(0);
    Theta_.min(100);

    if (debug)
    {
        Info<< "    max(Theta) = " << max(Theta_).value() << endl;
    }
}


// ************************************************************************* //
