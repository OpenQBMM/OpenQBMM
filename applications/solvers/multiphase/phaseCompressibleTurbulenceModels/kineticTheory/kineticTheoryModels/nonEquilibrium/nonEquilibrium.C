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
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(nonEquilibrium, 0);
    addToRunTimeSelectionTable
    (
        kineticTheoryModel,
        nonEquilibrium,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonEquilibrium::nonEquilibrium
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    kineticTheoryModel(dict, phase),
    alphaTheta_
    (
        "alphaTheta",
        dimless,
        dict
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::nonEquilibrium::~nonEquilibrium()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModels::nonEquilibrium::solve
(
    const volScalarField& beta,
    const volScalarField& alpha,
    const volTensorField& gradU,
    const volSymmTensorField D
)
{
    // Local references
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = phase_.alphaRhoPhi();
    const volVectorField& U = phase_.U();
    const volVectorField& Uc = phase_.fluid().otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    const volScalarField& da = phase_.d();

    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

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
          - fvc::ddt(alpha, rho, Theta_)
          + fvm::div
            (
                this->h2f()*alphaRhoPhi,
                Theta_,
                "div(" + alphaRhoPhi.name() + "," + Theta_.name() + ")"
            )
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
      + alpha*rho*particleTurbulenceModel.epsilon()

      + fvOptions(alpha, rho, Theta_)
    );

    ThetaEqn.relax();
    fvOptions.constrain(ThetaEqn);
    ThetaEqn.solve();
    fvOptions.correct(Theta_);

    Theta_.max(0);
    Theta_.min(100);

    if (debug)
    {
        Info<< "    max(Theta) = " << max(Theta_).value() << endl;
    }
}


// ************************************************************************* //
