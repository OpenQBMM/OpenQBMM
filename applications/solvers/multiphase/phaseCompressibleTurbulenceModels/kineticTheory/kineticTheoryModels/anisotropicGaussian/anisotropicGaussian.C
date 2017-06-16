/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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

#include "anisotropicGaussian.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "twoPhaseSystem.H"
#include "dragModel.H"
#include "fvm.H"
#include "fvc.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(anisotropicGaussian, 0);
    addToRunTimeSelectionTable
    (
        kineticTheoryModel,
        anisotropicGaussian,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::kineticTheoryModels::anisotropicGaussian::updateh2Fn()
{

    const dimensionedScalar smallPpk("small",dimensionSet(0, 2, -2, 0, 0), SMALL);

    g0_ = radialModel_->g0(phase_, alphaMinFriction_, alphaMax_);

    // This calculates the h2 function for the dense regime transport.
    if(h2FnMethod_.match("alphaG0"))
    {
        h2Fn_ = 1.0 - 1.0/(1.0 + sqr(phase_)*pow(g0_,h2FnParaPow_));
    }
    else if(h2FnMethod_.match("particlePressure"))
    {
        volScalarField ppk(max(phase_*Theta_,smallPpk));
        volScalarField pps(4.0*eta_*phase_*g0_*ppk + ppfr_);
        h2Fn_ = pow(pps/(pps + ppk),h2FnParaPow_);
    }
    else
    {
        FatalErrorIn("kineticTheoryModel::updateh2Fn: invalid h2FnMethod") << abort(FatalError);
    }

    h2Fn_.correctBoundaryConditions();

    return;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::anisotropicGaussian::anisotropicGaussian
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    kineticTheoryModel(dict, phase),
    alphaTheta_
    (
        "alphaTheta",
        dimensionSet(0, 0, 0, 0, 0),
        dict.lookup("alphaTheta")
    ),
    eta_(0.5*(1+e_)),
    ppfr_
    (
       IOobject
       (
           IOobject::groupName("ppfr", phase_.name()),
           phase.mesh().time().timeName(),
           phase.mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       phase.mesh(),
       dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),

    h2Fn_
    (
        IOobject
        (
            "h2Fn",
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       phase.mesh(),
       1.0
    ),

    Sigma_
    (
        IOobject
        (
            IOobject::groupName("Sigma", phase.name()),
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        2.0*nu_*dev(twoSymm(fvc::grad(phase_.U()))),
        Theta_.boundaryField().types()
    ),
    AGtransport_(phase.mesh(), dict, phase, Theta_, Sigma_),
    h2FnMethod_(dict.lookup("h2FnMethod")),
    h2FnParaPow_(readScalar(dict.lookup("h2FnParaPow", 2)))
{
    kappa_.dimensions().reset(kappa_.dimensions()/dimDensity);
    lambda_ =
        (8.0/3.0)/sqrt(constant::mathematical::pi)
       *phase_.d()*eta_*sqr(phase_)*g0_*sqrt(Theta_);

    ppfr_ =
        frictionalStressModel_->frictionalPressure
        (
            phase_,
            alphaMinFriction_,
            alphaMax_
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::anisotropicGaussian::~anisotropicGaussian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModels::anisotropicGaussian::updateViscosities()
{

    // Local references
    const dimensionedScalar smallRT("small",dimensionSet(0, 0, -1, 0, 0), SMALL);
    volScalarField alpha = max(phase_, scalar(0));
    const volScalarField& rho = phase_.rho();
    const volScalarField& da = phase_.d();
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    volScalarField alphaSqr(sqr(alpha));
    volScalarField thetaSqrt(sqrt(Theta_));
    volScalarField Kd
    (
        refCast<const twoPhaseSystem>(phase_.fluid()).drag(phase_).K()
    );

    g0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    // bulk viscosity
    lambda_ = (8.0/3.0)/sqrtPi*da*eta_*alphaSqr*g0_*thetaSqrt;

    volScalarField rTaupAlpha("rTaupAlpha", Kd/rho + smallRT );
    volScalarField rTaucAlpha
    (
        "rTaucAlpha",
        (6.0/sqrtPi/da)*g0_*max(alphaSqr, residualAlpha_)*thetaSqrt
    );

    // Particle viscosity
    nu_ =
        0.5*alphaSqr*(1.0 + (8.0/5.0)*eta_*(3*eta_ - 2.0)*alpha*g0_)
       *(h2Fn_ + (8.0/5.0)*eta_*alpha*g0_)
       *Theta_/(rTaupAlpha + eta_*(2.0-eta_)*rTaucAlpha)
      + (3.0/5.0)*lambda_;

    // Frictional pressure
    ppfr_ = frictionalStressModel_->frictionalPressure
    (
        phase_,
        alphaMinFriction_,
        alphaMax_
    );

    nuFric_ = frictionalStressModel_->nu
    (
        phase_,
        alphaMinFriction_,
        alphaMax_,
        ppfr_/rho,
        symm(fvc::grad(phase_.U()))
    );

    // Limit viscosity and add frictional viscosity
    nu_.min(maxNut_);
    nuFric_ = min(nuFric_, maxNut_ - nu_);
}


void Foam::kineticTheoryModels::anisotropicGaussian::correct()
{
    // Constants
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const dimensionedScalar smallRT
    (
        "small",
        dimensionSet(0, 0, -1, 0, 0),
        SMALL
    );

    // Refrences
    const volScalarField& alpha = phase_;
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaPhi = phase_.alphaPhi();
    const volScalarField& da = phase_.d();


    // Particle Strain-rate tensor
    tmp<volTensorField> tgradU(fvc::grad(phase_.U()));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));
    volSymmTensorField Sp(D - (1.0/3.0)*tr(D)*I);

    volScalarField alphaSqr(sqr(alpha));
    volScalarField thetaSqrt(sqrt(Theta_));

    g0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    // bulk viscosity
    lambda_ = (8.0/3.0)/sqrtPi*da*eta_*alphaSqr*g0_*thetaSqrt;

    // particle pressure, move bulk viscosity part to tau
    volScalarField  PsCoeff = alpha*(h2Fn_ + 4.0*eta_*alpha*g0_);
    volScalarField rTaupAlpha
    (
        "rTaupAlpha",
        phase_.fluid().Kd()/rho + smallRT
    );
    volScalarField rTaucAlphaCoeff
    (
        (6.0/sqrtPi/da)*g0_*max(alphaSqr, residualAlpha_)
    );
    volScalarField rTaucAlpha("rTaucAlpha", rTaucAlphaCoeff*thetaSqrt);

    // Particle viscosity
    volScalarField nutCoeff
    (
        0.5*alphaSqr*(1.0 + (8.0/5.0)*eta_*(3*eta_ - 2.0)*alpha*g0_)
       *h2Fn_*(1.0 + (8.0/5.0)*eta_*alpha*g0_)
    );

    nu_ =
        nutCoeff*Theta_/(rTaupAlpha + eta_*(2.0 - eta_)*rTaucAlpha)
      + (3.0/5.0)*lambda_;

    volSymmTensorField tau(2.0*nu_*Sp + lambda_*tr(D)*I);

    // 'thermal' conductivity
    volScalarField kappaCoeff
    (
        2.5*alphaSqr
       *(
            1.0
          + (12.0/5.0)*sqr(eta_)*(4.0*eta_ - 3.0)*alpha*g0_
        )
       *h2Fn_*(1.0 + (12.0/5.0)*eta_*alpha*g0_)
    );
    kappa_ =
        kappaCoeff*Theta_
       /(
            3.0*rTaupAlpha
          + 4.0*eta_*(41.0-33.0*eta_)*rTaucAlpha
        )
      + (3.0/2.0)*lambda_ ;

    // update collisional pressure
    volScalarField ppc(4.0*eta_*alphaSqr*g0_*Theta_ - lambda_*tr(D)) ;

    volSymmTensorField S2flux
    (
        (
            (h2Fn_*alpha*Theta_ + ppc)*Sp
          - nu_*(twoSymm(Sp & gradU) - (2.0/3.0)*(Sp && gradU)*I)
        )*2.0
    );
    volScalarField particleContinuityErr
    (
        fvc::ddt(alpha)
      + fvc::div(alphaPhi)
    );

    fvSymmTensorMatrix SigmaEqn
    (
        fvm::ddt(alpha, Sigma_)
      + fvm::div(alphaPhi, Sigma_)
      - fvc::Sp(particleContinuityErr, Sigma_)
      - fvm::laplacian(2.0/3.0*kappa_, Sigma_)
     ==
        S2flux
      + fvm::Sp
        (
          - (2.0*rTaupAlpha + (3.0 - e_)*(1.0 + e_)/2.0*rTaucAlpha),
            Sigma_
        )

    );

    SigmaEqn.relax();
    SigmaEqn.solve();

    // Construct the granular temperature equation
    fvScalarMatrix ThetaEqn
    (

        fvm::ddt(alpha, Theta_)
      + fvm::div(alphaPhi, Theta_)
      - fvc::Sp(particleContinuityErr, Theta_)
      - fvm::laplacian(2.0/3.0*kappa_, Theta_)
     ==
        fvm::SuSp(-2.0/3.0*((PsCoeff*I) && gradU), Theta_)
      + 2.0/3.0*(tau && gradU)
      + fvm::Sp
        (
          - (2.0*rTaupAlpha + (1.0 - sqr(e_))*rTaucAlpha),
            Theta_
        )
    );

    ThetaEqn.relax();
    ThetaEqn.solve();

    Theta_.max(0);
    Theta_.min(100);

    Theta_.correctBoundaryConditions();

    thetaSqrt = sqrt(Theta_);

    // update bulk viscosity
    lambda_ = (8.0/3.0)/sqrtPi*da*eta_*alphaSqr*g0_*thetaSqrt;

    // update Particle viscosity
    nu_ =
        nutCoeff*Theta_
       /(
            rTaupAlpha
          + eta_*(2.0-eta_)*rTaucAlphaCoeff*thetaSqrt
        )
      + (3.0/5.0)*lambda_;

    // Frictional pressure
    volScalarField ppfr = frictionalStressModel_->frictionalPressure
    (
        phase_,
        alphaMinFriction_,
        alphaMax_
    );

    // Update frictional shear viscosity
    nuFric_ = frictionalStressModel_->nu
    (
        phase_,
        alphaMinFriction_,
        alphaMax_,
        ppfr/rho,
        Sp
    );

    if (debug)
    {
        Info<< "    max(Theta) = " << max(Theta_).value() << nl
            << "    min(Theta) = " << min(Theta_).value() << nl
            << "    max(Sigma) = " << max(mag(Sigma_)).value() << nl
            << "    max(nut) = " << max(nu_).value() << endl;
    }
}

void Foam::kineticTheoryModels::anisotropicGaussian::transportMoments()
{
    updateh2Fn();
    surfaceScalarField h2Fnf = fvc::interpolate(h2Fn_);

    AGtransport_.solve(h2Fnf);

    surfaceScalarField& phi =
        phase_.mesh().lookupObjectRef<surfaceScalarField>
        (
            IOobject::groupName("phi", phase_.name())
        );
    phi = fvc::flux(phase_.U());
    correct_ = true;
}


Foam::scalar
Foam::kineticTheoryModels::anisotropicGaussian::maxUxDx() const
{
    return AGtransport_.maxUxDx();
}


// ************************************************************************* //
