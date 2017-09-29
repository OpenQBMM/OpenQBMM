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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class baseModel>
void Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::updateh2Fn()
{
    this->g0_ = this->radialModel_->g0
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
        this->frictionalStressModel_->frictionalPressure
        (
            this->phase_,
            this->alphaMinFriction_,
            this->alphaMax_
        ),
        this->e_
    );
    h2Fn_.max(this->residualAlpha_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class baseModel>
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::anisotropicGaussian
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    baseModel(dict, phase),
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
    eta_(0.5*(1.0 + this->e_)),
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
        2.0*this->nu_*dev(twoSymm(fvc::grad(this->phase_.U()))),
        this->Theta_.boundaryField().types()
    ),
    AGtransport_(phase.mesh(), dict, phase, this->Theta_, Sigma_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class baseModel>
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::
~anisotropicGaussian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class baseModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::Sigma() const
{
    return Sigma_;
}


template<class baseModel>
void Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::solve
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

    const scalar sqrtPi = sqrt(constant::mathematical::pi);

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
          + fvm::div
            (
                this->h2f()*alphaRhoPhi,
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

    baseModel::solve(beta, alpha, gradU, D);

    if (debug)
    {
        Info<< "    max(Sigma) = " << max(mag(Sigma_)).value() << endl;
    }
}


template<class baseModel>
void
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::transportMoments()
{
    Info<< "Transporting moments in dilute regime" << endl;

    updateh2Fn();
    AGtransport_.solve(this->h2f());
}


template<class baseModel>
Foam::scalar
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::maxUxDx() const
{
    return AGtransport_.maxUxDx();
}


template<class baseModel>
Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::h2() const
{
    return h2Fn_;
}


template<class baseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::h2f() const
{
    return fvc::interpolate(h2Fn_);
}


template<class baseModel>
Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::anisotropicGaussian<baseModel>::
ddtAlphaDilute() const
{
    return AGtransport_.ddtAlphaDilute();
}


// ************************************************************************* //
