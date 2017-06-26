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
    g0_ = radialModel_->g0(phase_, alphaMinFriction_, alphaMax_);
    h2Fn_ == h2Function_->h2
    (
        phase_,
        Theta_,
        g0_,
        phase_.rho(),
        phase_.d(),
        ppfr_,
        e_
    );
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
    alphaSigma_
    (
        "alphaSigma",
        dimensionSet(0, 0, 0, 0, 0),
        dict.lookup("alphaSigma")
    ),
    eta_(0.5*(1.0 + e_)),
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

    h2Function_(fluxSplittingFunction::New(dict)),
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
    AGtransport_(phase.mesh(), dict, phase, Theta_, Sigma_)
{
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

void Foam::kineticTheoryModels::anisotropicGaussian::correct()
{
    // Local references
    volScalarField alpha(max(phase_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = phase_.alphaRhoPhi();
    const volVectorField& U = phase_.U();
    const volVectorField& Uc = phase_.fluid().otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    const volScalarField& da = phase_.d();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));
    volSymmTensorField Sp(D - (1.0/3.0)*tr(D)*I);

    // Calculating the radial distribution function
    g0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));
    volScalarField alphaSqr("alphaSqr", sqr(alpha));

    // Bulk viscosity  p. 45 (KongFox et al. 1984).
    lambda_ = (4.0/3.0)*alphaSqr*da*g0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    // Dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
    (
        "gammaCoeff",
        12.0*(1.0 - sqr(e_))
        *max(alphaSqr, residualAlpha_)
        *rho*g0_*(1.0/da)*ThetaSqrt/sqrtPi
    );

    // Drag
    volScalarField beta
    (
        refCast<const twoPhaseSystem>(phase_.fluid()).drag(phase_).K()
    );
    volScalarField rTauc
    (
        "rTauc",
        6.0*sqrt(Theta_)*max(alpha, phase_.residualAlpha())*g0_/(da*sqrtPi)
    );

    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1("J1", 3.0*beta);
    volScalarField J2
    (
        "J2",
        0.25*sqr(beta)*da*magSqr(U - Uc)
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
            g0_,
            rho,
            e_
        )
    );

    // Particle viscosity
    nu_ = (1.0 + 8.0/5.0*eta_*alpha*g0_)*max(h2Fn_, residualAlpha_)
       *viscosityModel_->nu(phase_, Theta_, g0_, rho, da, e_)
      + 0.6*lambda_;

    // 'thermal' conductivity
    kappa_ = (1.0 + 12.0/5.0*eta_*alpha*g0_)*max(h2Fn_, residualAlpha_)
       *conductivityModel_->kappa(phase_, Theta_, g0_, rho, da, e_)
      + 1.5*lambda_*rho;

    // Stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau
    (
        rho*(2.0*nu_*D + (lambda_ - (2.0/3.0)*nu_)*tr(D)*I)
    );

    fv::options& fvOptions(fv::options::New(phase_.fluid().mesh()));
    const PhaseCompressibleTurbulenceModel<phaseModel>&
        particleTurbulenceModel =
            U.db().lookupObject<PhaseCompressibleTurbulenceModel<phaseModel> >
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    phase_.name()
                )
            );

    // Solve Sigma equation (2nd order moments)
    {
        nu_ += nuFric_;

        volSymmTensorField S2flux
        (
            "S2flux",
            2.0
           *(
                h2Fn_*alpha*rho*Theta_
              + 4.0*rho*eta_*alphaSqr*g0_*Theta_
              - rho*lambda_*tr(D)
            )*Sp
          - 2.0*rho*alpha*nu_
           *(
               twoSymm(Sp & gradU)
             - (2.0/3.0)*(Sp && gradU)*I
            )
        );

        fvSymmTensorMatrix SigmaEqn
        (
            fvm::ddt(alpha, rho, Sigma_)
          - fvc::ddt(alpha, rho, Sigma_)
          + fvm::div
            (
                hydrodynamicScalef(alphaRhoPhi),
                Sigma_,
                "div(" + alphaRhoPhi.name() + "," + Sigma_.name() + ")"
            )
          - fvc::Sp
            (
                fvc::ddt(alpha, rho)
              + fvc::div(alphaRhoPhi),
                Sigma_
            )
          - fvm::laplacian
            (
                kappa_ + rho*particleTurbulenceModel.nut()/alphaSigma_,
                Sigma_,
                "laplacian(kappa,Sigma)"
            )
         ==
            S2flux
          - fvm::Sp
            (
                alpha*(2.0*beta + (3.0 - e_)*(1.0 + e_)/2.0*rTauc*rho),
                Sigma_
            )
        );

        SigmaEqn.relax();
        SigmaEqn.solve();
    }

    // Construct the granular temperature equation (Eq. 3.20, p. 44)
    // NB. note that there are two typos in Eq. 3.20:
    //     Ps should be without grad
    //     the laplacian has the wrong sign
    fvScalarMatrix ThetaEqn
    (
        1.5*
        (
            fvm::ddt(alpha, rho, Theta_)
          - fvc::ddt(alpha, rho, Theta_)
          + fvm::div(alphaRhoPhi, Theta_)
          - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Theta_)
        )
      - fvm::laplacian
        (
            kappa_ + rho*particleTurbulenceModel.nut()/alphaTheta_,
            Theta_,
            "laplacian(kappa,Theta)"
        )
     ==
      - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
      + (tau && gradU)
      + fvm::Sp(-gammaCoeff, Theta_)
      + fvm::Sp(-J1, Theta_)
      + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
      + alpha*rho*particleTurbulenceModel.epsilon()

      + fvOptions(alpha, rho, Theta_)
    );

    ThetaEqn.relax();
    fvOptions.constrain(ThetaEqn);
    ThetaEqn.solve();
    fvOptions.correct(Theta_);

    Theta_.max(0);
    Theta_.min(100);

    ThetaSqrt = sqrt(Theta_);

    {
        // Bulk viscosity  p. 45 (KongFox et al. 1984).
        lambda_ = (4.0/3.0)*alphaSqr*da*g0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

        // particle viscosity (Table 3.2, p.47)
        nu_ = viscosityModel_->nu(phase_, Theta_, g0_, rho, da, e_);

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
            D
        );

        // Limit viscosity and add frictional viscosity
        nu_.min(maxNut_);
        nuFric_ = min(nuFric_, maxNut_ - nu_);
    }

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
    Info<< "Transporting moments" << endl;
    updateh2Fn();
    surfaceScalarField h2Fnf = fvc::interpolate(h2Fn_);

    AGtransport_.solve(h2Fnf);

    surfaceScalarField& phi =
        phase_.mesh().lookupObjectRef<surfaceScalarField>(phase_.phi().name());
    phi = fvc::flux(phase_.U());
}


Foam::scalar
Foam::kineticTheoryModels::anisotropicGaussian::maxUxDx() const
{
    return AGtransport_.maxUxDx();
}


// ************************************************************************* //
