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
    if (!correct_)
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
    else
    {
        // Local references
        volScalarField alpha = max(phase_, scalar(0));
        const volScalarField& rho = phase_.rho();
        const surfaceScalarField& alphaRhoPhi = phase_.alphaRhoPhi();
        const volVectorField& U = phase_.U();
        const volVectorField& Uc = phase_.fluid().otherPhase(phase_).U();
        const volScalarField& da = phase_.d();

        volScalarField alphaSqr(sqr(alpha));
        volScalarField thetaSqrt(sqrt(Theta_));

        // Particle Strain-rate tensor
        tmp<volTensorField> tgradU(fvc::grad(phase_.U()));
        const volTensorField& gradU(tgradU());
        volSymmTensorField D(symm(gradU));
        volSymmTensorField Sp(D - (1.0/3.0)*tr(D)*I);

        // Drag
        volScalarField beta
        (
            refCast<const twoPhaseSystem>(phase_.fluid()).drag(phase_).K()
        );

        // Constants
        const scalar sqrtPi = sqrt(constant::mathematical::pi);
        dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
        dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));
        const dimensionedScalar smallRT
        (
            "small",
            dimensionSet(0, 0, -1, 0, 0),
            SMALL
        );

        // Calculating the radial distribution function
        g0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

        volScalarField rTaupAlpha("rTaupAlpha", beta/rho + smallRT );

        volScalarField rTaucAlphaCoeff
        (
            (6.0/sqrtPi/da)*g0_*max(alphaSqr, phase_.residualAlpha())
        );

        volScalarField rTaucAlpha("rTaucAlpha", rTaucAlphaCoeff*thetaSqrt);

        // Particle viscosity
        volScalarField nutCoeff
        (
            0.5*alphaSqr*(1.0 + (8.0/5.0)*eta_*(3*eta_ - 2.0)*alpha*g0_)
            *h2Fn_*( 1.0 + (8.0/5.0)*eta_*alpha*g0_)
        );

        // 'thermal' conductivity
        volScalarField kappaCoeff
        (
            2.5*alphaSqr*(1.0 + (12.0/5.0)*sqr(eta_)
           *(4.0*eta_ - 3.0)*alpha*g0_)
           *h2Fn_*(1.0 + (12.0/5.0)*eta_*alpha*g0_)
        );
        kappa_ =
            kappaCoeff*Theta_
           /(3.0*rTaupAlpha + 4.0*eta_*(41.0-33.0*eta_)*rTaucAlpha)
          + (3.0/2.0)*lambda_;


        // update collisional pressure
        volScalarField ppc(4.0*eta_*alphaSqr*g0_*Theta_ - lambda_*tr(D));

        nu_ =
            nutCoeff*Theta_/(rTaupAlpha + eta_*(2.0-eta_)*rTaucAlpha)
          + (3.0/5.0)*lambda_;

        // Stress tensor, Definitions, Table 3.1, p. 43
        volSymmTensorField tau(rho*(2.0*nu_*Sp + lambda_*tr(D)*I));

        volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*g0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

//         // Stress tensor, Definitions, Table 3.1, p. 43
//         volSymmTensorField tau
//         (
//             rho*(2.0*nu_*D + (lambda_ - (2.0/3.0)*nu_)*tr(D)*I)
//         );

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

        // particle pressure, move bulk viscosity part to tau
        volScalarField  PsCoeff = alpha*(h2Fn_ + 4.0*eta_*alpha*g0_);
//         // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
//         volScalarField PsCoeff
//         (
//             granularPressureModel_->granularPressureCoeff
//             (
//                 alpha,
//                 g0_,
//                 rho,
//                 e_
//             )
//         );

        // 'thermal' conductivity (Table 3.3, p. 49)
        kappa_ = conductivityModel_->kappa(alpha, Theta_, g0_, rho, da, e_);

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

        volSymmTensorField S2flux
        (
            (
                ppc*Sp
              - nu_*alpha*rho
               *(
                   twoSymm(Sp & gradU)
                 - (2.0/3.0)*(Sp && gradU)*I
                )
            )*2.0
        );

        fvSymmTensorMatrix SigmaEqn
        (
            fvm::ddt(alpha, rho, Sigma_)
          + fvm::div(alphaRhoPhi*fvc::interpolate(h2Fn_), Sigma_)
          - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Sigma_)
          - fvm::laplacian(2.0/3.0*kappa_, Sigma_)
         ==
            S2flux
          - fvm::Sp
            (
                (
                    2.0*rTaupAlpha
                  + (3.0-e_)*(1.0+e_)/2.0*rTaucAlpha
                )*alpha*rho,
                Sigma_
            )

        );

        SigmaEqn.relax();
        SigmaEqn.solve();

        if (debug)
        {
            Info<< typeName << ':' << nl
                << "    max(Theta) = " << max(Theta_).value() << nl
                << "    max(Sigma) = " << max(mag(Sigma_)).value() << endl;
        }
    }
}




// ************************************************************************* //
