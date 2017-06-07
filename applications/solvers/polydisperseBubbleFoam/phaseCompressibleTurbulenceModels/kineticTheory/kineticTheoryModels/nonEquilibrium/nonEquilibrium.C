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

#include "nonEquilibrium.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "twoPhaseSystem.H"
#include "dragModel.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"


#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
defineTypeNameAndDebug(nonEquilibrium, 0);
addToRunTimeSelectionTable(nonEquilibrium, nonEquilibrium, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonEquilibrium::nonEquilibrium
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    nonEquilibrium(dict, phase),
    alphaTheta_
    (
        "alphaTheta",
        dimless,
        dict.lookup("alphaTheta")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonEquilibrium::~nonEquilibrium()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModels::nonEquilibrium::correct()
{
    // Local references
    volScalarField alpha(max(alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;
    const volVectorField& Uc_ =
        refCast<const twoPhaseSystem>(phase_.fluid()).otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    // Particle viscosity (Table 3.2, p.47)
    nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

    // Bulk viscosity  p. 45 (Lun et al. 1984).
    lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    // Stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau
    (
        rho*(2.0*nut_*D + (lambda_ - (2.0/3.0)*nut_)*tr(D)*I)
    );

    // Dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
    (
        "gammaCoeff",
        12.0*(1.0 - sqr(e_))
        *max(sqr(alpha), residualAlpha_)
        *rho*gs0_*(1.0/da)*ThetaSqrt/sqrtPi
    );

    // Drag
    volScalarField beta
    (
        refCast<const twoPhaseSystem>(phase_.fluid()).drag(phase_).K()
    );

    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1("J1", 3.0*beta);
    volScalarField J2
    (
        "J2",
        0.25*sqr(beta)*da*magSqr(U - Uc_)
        /(
            max(alpha, residualAlpha_)*rho
            *sqrtPi*(ThetaSqrt + ThetaSmallSqrt)
        )
    );

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff
    (
        granularPressureModel_->granularPressureCoeff
        (
            alpha,
            gs0_,
            rho,
            e_
        )
    );

    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);

    fv::options& fvOptions(fv::options::New(mesh_));

    // Construct the granular temperature equation (Eq. 3.20, p. 44)
    // NB. note that there are two typos in Eq. 3.20:
    //     Ps should be without grad
    //     the laplacian has the wrong sign
    fvScalarMatrix ThetaEqn
    (
        1.5*
        (
            fvm::ddt(alpha, rho, Theta_)
            + fvm::div(alphaRhoPhi, Theta_)
            - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Theta_)
        )
        - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
        ==
        - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
        + (tau && gradU)
        + fvm::Sp(-gammaCoeff, Theta_)
        + fvm::Sp(-J1, Theta_)
        + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
        + fvOptions(alpha, rho, Theta_)
    );

    ThetaEqn.relax();
    fvOptions.constrain(ThetaEqn);
    ThetaEqn.solve();
    fvOptions.correct(Theta_);

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }
}


// ************************************************************************* //
