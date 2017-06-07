/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "kineticTheoryModel.H"
#include "mathematicalConstants.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::kineticTheoryModel
(
    const dictionary dict,
    const phaseModel& phase
)
:
    kineticTheoryModel(dict, phase),
    mesh_(mesh),
    particles_(particles),
    gas_(gas),

    eta_(0.5*(1+e_)),
    alphap_(particles_),

    kappa_
    (
       IOobject
       (
           IOobject::groupName("kappa", particles_.name()),
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh_,
       dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),
    ppfr_
    (
       IOobject
       (
           IOobject::groupName("ppfr", particles_.name()),
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh_,
       dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
    ),

    h2Fn_
    (
        IOobject
        (
            "h2Fn",
            runTime.timeName(),
            mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       mesh,
       1.0
    ),

    Sigma_
    (
       IOobject
       (
           IOobject::groupName("Sigma", particles_.name()),
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::AUTO_WRITE
       ),
       mesh_
    ),
    diluteAGmodel_(mesh, *this, particles, Theta_, Sigma_),
    h2FnMethod_(this->lookup("h2FnMethod")),
    h2FnParaPow_(readLabel(this->lookup("h2FnParaPow")))

{

    g0_ = radialModel_->g0(alphap_, alphaMinFriction_, alphaMax_);

    lambda_ =
        (8.0/3.0)/sqrt(constant::mathematical::pi)
       *particles_.d()*eta_*sqr(alphap_)*gs0_*sqrt(Theta_);

    ppfr_ = frictionalStressModel_->frictionalPressure(alphap_);


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModel::pPrime() const
{

    return  (h2Fn_ + 4.0*alphap_*eta_*(2.0*radialModel_->g0(alphap_, alphaMinFriction_, alphaMax_)
			+ radialModel_->g0prime(alphap_, alphaMinFriction_, alphaMax_)*alphap_))*Theta_
            + frictionalStressModel_->frictionalPressurePrime(alphap_);

}


void Foam::kineticTheoryModel::updateh2Fn()
{

    const dimensionedScalar smallPpk("small",dimensionSet(0, 2, -2, 0, 0), SMALL);

    gs0_ = radialModel_->g0(alphap_, alphaMinFriction_, alphaMax_);

    // This calculates the h2 function for the dense regime transport.
    if(h2FnMethod_.match("alphaG0"))
    {
        h2Fn_ = 1.0 - 1.0/(1.0 + sqr(alphap_)*pow(gs0_,h2FnParaPow_));
    }
    else if(h2FnMethod_.match("particlePressure"))
    {
        volScalarField ppk(max(alphap_*Theta_,smallPpk));
        volScalarField pps(4.0*eta_*alphap_*gs0_*ppk + ppfr_);
        h2Fn_ = pow(pps/(pps + ppk),h2FnParaPow_);
    }
    else
    {
        FatalErrorIn("kineticTheoryModel::updateh2Fn: invalid h2FnMethod") << abort(FatalError);
    }

    h2Fn_.correctBoundaryConditions();

    return;
}


void Foam::kineticTheoryModel::solveDilute(const surfaceScalarField& h2f)
{
    if(diluteAGmodel_.AGmodel() && (h2FnParaPow_>0))  diluteAGmodel_.solve(h2f);
}


void Foam::kineticTheoryModel::updateViscosity()
{

    // Local references
    const dimensionedScalar smallRT("small",dimensionSet(0, 0, -1, 0, 0), SMALL);
    const dimensionedScalar& rho = phase_.rho();
    const dimensionedScalar& da = phase_.d();
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    volScalarField alphaSqr(sqr(alphap_));
    volScalarField thetaSqrt(sqrt(Theta_));
    volScalarField Kd = phase_.fluid().drag(phase_).K();

    g0_ = radialModel_->g0(alphap_, alphaMinFriction_, alphaMax_);

    // bulk viscosity
    lambda_ = (8.0/3.0)/sqrtPi*da*eta_*alphaSqr*gs0_*thetaSqrt;

    volScalarField rTaupAlpha("rTaupAlpha", Kd/rho + smallRT );
    volScalarField rTaucAlpha("rTaucAlpha", (6.0/sqrtPi/da)*gs0_*max(alphaSqr, residualAlpha_)*thetaSqrt);

    // Particle viscosity
    nu_ =  0.5*alphaSqr*(1.0 + (8.0/5.0)*eta_*(3*eta_ - 2.0)*alphap_*gs0_)
            *(h2Fn_ + (8.0/5.0)*eta_*alphap_*gs0_)*Theta_/(rTaupAlpha + eta_*(2.0-eta_)*rTaucAlpha) + (3.0/5.0)*lambda_;

    // Frictional pressure
    ppfr_ = frictionalStressModel_->frictionalPressure(alphap_);

    tmp<volTensorField> tgradU(fvc::grad(particles_.U()));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));  // D = 0.5*(gradU + gradU^T)
    volSymmTensorField Sp(D - (1.0/3.0)*tr(D)*I);
}


void Foam::kineticTheoryModel::solveDense
(
    const volScalarField& K,
    const volScalarField& particleContinuityErr
)
{

    // Local references
    const dimensionedScalar smallRT("small",dimensionSet(0, 0, -1, 0, 0), SMALL);
    const dimensionedScalar& rho = particles_.rho();
    const surfaceScalarField& alphaPhi = particles_.alphaPhi();
    const dimensionedScalar& da = particles_.d();
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    // Particle Strain-rate tensor
    tmp<volTensorField> tgradU(fvc::grad(particles_.U()));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));  // D = 0.5*(gradU + gradU^T)
    volSymmTensorField Sp(D - (1.0/3.0)*tr(D)*I);

    gs0_ = radialModel_->g0(alphap_, alphaMinFriction_, alphaMax_);

    volScalarField alphaSqr(sqr(alphap_));

    volScalarField thetaSqrt(sqrt(Theta_));

    // bulk viscosity
    lambda_ = (8.0/3.0)/sqrtPi*da*eta_*alphaSqr*gs0_*thetaSqrt;

    // particle pressure, move bulk viscosity part to tau
    volScalarField  PsCoeff = alphap_*(h2Fn_ + 4.0*eta_*alphap_*gs0_);

    volScalarField rTaupAlpha("rTaupAlpha", K/rho + smallRT );

    volScalarField rTaucAlphaCoeff((6.0/sqrtPi/da)*gs0_*max(alphaSqr, residualAlpha_));
    volScalarField rTaucAlpha("rTaucAlpha", rTaucAlphaCoeff*thetaSqrt);

    // Particle viscosity
    volScalarField nutCoeff
    (
        0.5*alphaSqr*(1.0 + (8.0/5.0)*eta_*(3*eta_ - 2.0)*alphap_*gs0_)
        *h2Fn_*( 1.0 + (8.0/5.0)*eta_*alphap_*gs0_)
    );

    nut_ = nutCoeff*Theta_/(rTaupAlpha + eta_*(2.0-eta_)*rTaucAlpha) + (3.0/5.0)*lambda_;

    volSymmTensorField tau(2.0*nut_*Sp + lambda_*tr(D)*I);

    // 'thermal' conductivity
    volScalarField kappaCoeff
    (
        2.5*alphaSqr*(1.0 + (12.0/5.0)*sqr(eta_)*(4.0*eta_ - 3.0)*alphap_*gs0_)
        *h2Fn_*(1.0 + (12.0/5.0)*eta_*alphap_*gs0_)
    );
    kappa_ = kappaCoeff*Theta_/(3.0*rTaupAlpha + 4.0*eta_*(41.0-33.0*eta_)*rTaucAlpha) + (3.0/2.0)*lambda_ ;

    if (diluteAGmodel_.AGmodel() )
    {
        // update collisional pressure
        volScalarField ppc(4.0*eta_*alphaSqr*gs0_*Theta_ - lambda_*tr(D)) ;

        volSymmTensorField S2flux(2.0*((h2Fn_*alphap_*Theta_ + ppc)*Sp - nut_*( twoSymm(Sp & gradU) - (2.0/3.0)*(Sp && gradU)*I)));

        fvSymmTensorMatrix SigmaEqn
        (
            fvm::ddt(alphap_, Sigma_)
            + fvm::div(alphaPhi, Sigma_)
            - fvc::Sp(particleContinuityErr, Sigma_)

            - fvm::laplacian(2.0/3.0*kappa_, Sigma_)
            ==
            S2flux
            + fvm::Sp(-(2.0*rTaupAlpha + (3.0-e_)*(1.0+e_)/2.0*rTaucAlpha),Sigma_)

        );

        SigmaEqn.relax();
        SigmaEqn.solve();

        Sigma_.correctBoundaryConditions();

    }

    // Construct the granular temperature equation
    fvScalarMatrix ThetaEqn
    (

        fvm::ddt(alphap_, Theta_)
        + fvm::div(alphaPhi, Theta_)
        - fvc::Sp(particleContinuityErr, Theta_)

        - fvm::laplacian(2.0/3.0*kappa_, Theta_)

        ==

        fvm::SuSp(- 2.0/3.0*((PsCoeff*I) && gradU), Theta_)
        + 2.0/3.0*( tau && gradU)
        + fvm::Sp(- (2.0*rTaupAlpha + (1.0 - sqr(e_))*rTaucAlpha), Theta_)
    );

    ThetaEqn.relax();
    ThetaEqn.solve();

    Theta_.max(0);
    Theta_.min(100);

    Theta_.correctBoundaryConditions();

    thetaSqrt = sqrt(Theta_);

    // update bulk viscosity
    lambda_ = (8.0/3.0)/sqrtPi*da*eta_*alphaSqr*gs0_*thetaSqrt;

    // update Particle viscosity
    nut_ = nutCoeff*Theta_/(rTaupAlpha + eta_*(2.0-eta_)*rTaucAlphaCoeff*thetaSqrt) + (3.0/5.0)*lambda_;

    // Frictional pressure
    ppfr_ = frictionalStressModel_->frictionalPressure(alphap_);

    // Add frictional shear viscosity
    nut_ += frictionalStressModel_->nu(alphap_,ppfr_,Sp);

    // Limit viscosity
    nut_.min(100);

    if (debug)
    {
        Info<< "    max(Theta) = " << max(Theta_).value() << nl
            << "    min(Theta) = " << min(Theta_).value() << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }

    return ;

}




// ************************************************************************* //
