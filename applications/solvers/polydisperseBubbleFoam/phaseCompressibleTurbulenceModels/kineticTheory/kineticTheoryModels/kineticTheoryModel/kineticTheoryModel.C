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

#include "kineticTheoryModel.H"
#include "twoPhaseSystem.H"
#include "mathematicalConstants.H"
#include "fvOptions.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kineticTheoryModel, 0);
    defineRunTimeSelectionTable(kineticTheoryModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::kineticTheoryModel
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, phase.name()),
            phase.db().time().timeName(),
            phase.mesh()
        )
    ),

    phase_(phase),

    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            dict
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            dict
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            dict
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            dict
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            dict
        )
    ),

    equilibrium_(dict.lookup("equilibrium")),
    e_("e", dimless, dict),
    alphaMax_("alphaMax", dimless, dict),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        dict
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict
    ),

    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        dict.lookupOrDefault<scalar>("maxNut", 1000)
    ),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    g0_
    (
        IOobject
        (
            IOobject::groupName("g0", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),

    nu_
    (
        IOobject
        (
            IOobject::groupName("nu", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModel::nuEff() const
{
    return nu_ + nuFric_;
}


Foam::tmp<Foam::scalarField>
Foam::kineticTheoryModel::nuEff(const label patchi) const
{
    return nu_.boundaryField()[patchi] + nuFric_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModel::g0Prime() const
{
    return radialModel_->g0prime(phase_, alphaMinFriction_, alphaMax_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();

    tmp<volScalarField> tpPrime
    (
        Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            phase_,
            radialModel_->g0(phase_, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(phase_, alphaMinFriction_, alphaMax_),
            rho,
            e_
        )
     +  frictionalStressModel_->frictionalPressurePrime
        (
            phase_,
            alphaMinFriction_,
            alphaMax_
        )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::kineticTheoryModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


void Foam::kineticTheoryModel::correct()
{
    // Local references
    volScalarField alpha(max(phase_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = phase_.alphaRhoPhi();
    const volVectorField& U = phase_.U();
    const volVectorField& Uc_ =
        refCast<const twoPhaseSystem>(phase_.fluid()).otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    g0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    if (!equilibrium_)
    {
        // Particle viscosity (Table 3.2, p.47)
        nu_ = viscosityModel_->nu(alpha, Theta_, g0_, rho, da, e_);

        volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*g0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

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
                g0_,
                rho,
                e_
            )
        );

        // 'thermal' conductivity (Table 3.3, p. 49)
        kappa_ = conductivityModel_->kappa(alpha, Theta_, g0_, rho, da, e_);

        fv::options& fvOptions(fv::options::New(phase_.mesh()));

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
    }
    else
    {
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
           *fvc::div(phase_.phi())
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

        Theta_ = sqr
        (
            (l1 + sqrt(l2 + l3))
           /(2.0*max(alpha, residualAlpha_)*K4)
        );

        kappa_ = conductivityModel_->kappa(alpha, Theta_, g0_, rho, da, e_);
    }

    Theta_.max(0);
    Theta_.min(100);

    {
        // particle viscosity (Table 3.2, p.47)
        nu_ = viscosityModel_->nu(alpha, Theta_, g0_, rho, da, e_);

        volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*g0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

        // Frictional pressure
        volScalarField pf
        (
            frictionalStressModel_->frictionalPressure
            (
                phase_,
                alphaMinFriction_,
                alphaMax_
            )
        );

        nuFric_ = frictionalStressModel_->nu
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            pf/rho,
            D
        );

        // Limit viscosity
        nu_.min(maxNut_);
        nuFric_ = min(nuFric_, maxNut_ - nu_);
    }

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nu) = " << max(nu_).value() << endl;
    }
}

// ************************************************************************* //
