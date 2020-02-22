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
    this->g0_ = this->radialModel_->g0
    (
        this->phase_,
        this->alphaMinFriction_,
        this->alphaMax_
    );
    PsFric_ = this->frictionalStressModel_->frictionalPressure
    (
        this->phase_,
        this->alphaMinFriction_,
        this->alphaMax_
    );

    h2Fn_ = h2Function_->h2
    (
        this->phase_,
        this->Theta_,
        this->g0_,
        this->phase_.rho(),
        this->phase_.d(),
        PsFric_,
        this->e_
    );
    h2Fn_.max(this->residualAlpha_);
    h2Fn_.correctBoundaryConditions();
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
        dict
    ),
    alphaSigma_
    (
        "alphaSigma",
        dimensionSet(0, 0, 0, 0, 0),
        dict
    ),
    eta_(0.5*(1.0 + this->e_)),
    h2Function_(fluxSplittingFunction::New(dict)),
    h2Fn_
    (
        IOobject
        (
            IOobject::groupName("h2Fn", phase.name()),
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       phase.mesh(),
       1.0
    ),
    PsFric_
    (
        IOobject
        (
            IOobject::groupName("PsFric", phase.name()),
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       phase.mesh(),
       dimensionedScalar("zero", dimPressure, 0.0)
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
        phase.mesh()
    ),
    AGtransport_(phase.mesh(), dict, phase, this->Theta_, Sigma_)
{
    Sigma_ = 2.0*this->nu_*dev(twoSymm(fvc::grad(this->phase_.U())));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::anisotropicGaussian::
~anisotropicGaussian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::anisotropicGaussian::Sigma() const
{
    return Sigma_;
}


void Foam::kineticTheoryModels::anisotropicGaussian::solve
(
    const volScalarField& beta,
    const volScalarField& alpha,
    const volTensorField& gradU,
    const volSymmTensorField D
)
{
    // Local references
    const volScalarField& rho = this->phase_.rho();
    const surfaceScalarField& alphaRhoPhi = this->phase_.alphaRhoPhi();
    const volVectorField& U = this->phase_.U();
    const volVectorField& Uc = phase_.fluid().otherPhase(phase_).U();

    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

    const volScalarField& da = this->phase_.d();
    volSymmTensorField Sp(D - (1.0/3.0)*tr(D)*I);


    // Drag
    volScalarField rTauc
    (
        "rTauc",
        6.0*sqrt(this->Theta_)
       *max
        (
            alpha,
            this->phase_.residualAlpha()
        )*this->g0_/(da*sqrtPi)
    );

    // Particle viscosity
    this->nu_ *= (1.0 + 8.0/5.0*this->eta_*alpha*this->g0_)*h2Fn_;
    this->nu_ += 3.0/5.0*this->lambda_;

    // 'thermal' conductivity
    this->kappa_ *= (1.0 + 12.0/5.0*this->eta_*alpha*this->g0_)*h2Fn_;
    this->kappa_ += 3.0/2.0*this->lambda_*rho;

    fv::options& fvOptions(fv::options::New(this->phase_.fluid().mesh()));
//     const PhaseCompressibleTurbulenceModel<phaseModel>&
//         particleTurbulenceModel =
//             U.db().lookupObject<PhaseCompressibleTurbulenceModel<phaseModel> >
//             (
//                 IOobject::groupName
//                 (
//                     turbulenceModel::propertiesName,
//                     this->phase_.name()
//                 )
//             );

    // Solve Sigma equation (2nd order moments)
    {
        volSymmTensorField S2flux
        (
            "S2flux",
            2.0*Sp
           *(
                this->granularPressureModel_->granularPressureCoeff
                (
                    alpha,
                    this->g0_,
                    rho,
                    this->e_
                )*this->Theta_
              - rho*alpha*this->lambda_*tr(D)
            )
          - 2.0*rho*alpha*this->nu_
           *(
               twoSymm(Sp & gradU)
             - (2.0/3.0)*(Sp && gradU)*I
            )
        );

        fvSymmTensorMatrix SigmaEqn
        (
            fvm::ddt(alpha, rho, Sigma_)
          - fvc::ddt(alpha, rho, Sigma_)
          + fvm::div(alphaRhoPhi, Sigma_)
          - fvc::Sp
            (
                fvc::ddt(alpha, rho)
              + fvc::div(alphaRhoPhi),
                Sigma_
            )
          - fvm::laplacian
            (
                2.0/3.0*this->kappa_,
//               + rho*particleTurbulenceModel.nut()/alphaSigma_,
                Sigma_,
                "laplacian(kappa,Sigma)"
            )
         ==
            S2flux
          - fvm::Sp
            (
                alpha
               *(
                   2.0*beta + (3.0 - this->e_)*(1.0 + this->e_)/2.0*rTauc*rho
                ),
                Sigma_
            )
          + fvOptions(alpha, rho, Sigma_)
//           + alpha*rho*particleTurbulenceModel.epsilon()
        );

        SigmaEqn.relax();
        fvOptions.constrain(SigmaEqn);
        SigmaEqn.solve();
        fvOptions.correct(Sigma_);
    }

    // Solve granular temperature transport
    // Stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau
    (
        rho*(2.0*nu_*D + (lambda_ - (2.0/3.0)*nu_)*tr(D)*I)
    );

    // Dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
    (
        "gammaCoeff",
        12.0*(1.0 - sqr(e_))
        *max(sqr(alpha), residualAlpha_)
        *rho*g0_*(1.0/da)*ThetaSqrt/sqrtPi
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
            kappa_,
//           + rho*particleTurbulenceModel.nut()/alphaTheta_,
            Theta_,
            "laplacian(kappa,Theta)"
        )
     ==
      - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
      + (tau && gradU)
      + fvm::Sp(-gammaCoeff, Theta_)
      + fvm::Sp(-J1, Theta_)
      + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
//       + alpha*rho*particleTurbulenceModel.epsilon()

      + fvOptions(alpha, rho, Theta_)
    );

    ThetaEqn.relax();
    fvOptions.constrain(ThetaEqn);
    ThetaEqn.solve();
    fvOptions.correct(Theta_);

    Theta_.max(0);
    Theta_.min(100);

    PsFric_ =
        frictionalStressModel_->frictionalPressurePrime
        (
            phase_,
            alphaMinFriction_,
            alphaMax_
        );

    if (debug)
    {
        Info<< "    max(Sigma) = " << max(mag(Sigma_)).value()
            << "    max(Theta) = " << max(Theta_).value() << endl;
    }
}


void Foam::kineticTheoryModels::anisotropicGaussian::transportMoments()
{
    Info<< "Transporting moments in dilute regime" << endl;

    updateh2Fn();
    AGtransport_.solve(this->h2f());
}


Foam::scalar Foam::kineticTheoryModels::anisotropicGaussian::maxUxDx() const
{
    return AGtransport_.maxUxDx();
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::anisotropicGaussian::h2() const
{
    return h2Fn_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::kineticTheoryModels::anisotropicGaussian::h2f() const
{
    return fvc::interpolate(h2Fn_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::anisotropicGaussian::ddtAlphaDilute() const
{
    return AGtransport_.ddtAlphaDilute();
}


// ************************************************************************* //
