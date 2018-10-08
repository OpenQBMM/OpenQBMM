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
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
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
        dimensionSet(0, 2, -1, 0, 0),
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
            IOobject::groupName("gs0", phase.name()),
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
{
    g0_ = radialModel_->g0(phase_, alphaMinFriction_, alphaMax_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField> Foam::kineticTheoryModel::Sigma() const
{
    tmp<volTensorField> gradU(fvc::grad(phase_.U()));
    return nu_*(twoSymm(gradU()) - (2.0/3.0)*tr(gradU())*I);
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModel::nuEff() const
{
    tmp<volScalarField> nuEff(nu_ + nuFric_);
    nuEff.ref().min(100);
    return nuEff;
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


void Foam::kineticTheoryModel::update()
{
    g0_ = radialModel_->g0
    (
        max(phase_, scalar(0)),
        alphaMinFriction_,
        phase_.alphaMax()
    );

    // particle viscosity (Table 3.2, p.47)
    nu_ = viscosityModel_->nu
    (
        phase_,
        Theta_,
        g0_,
        phase_.rho(),
        phase_.d(),
        e_
    );

    // Bulk viscosity  p. 45 (Lun et al. 1984).
    lambda_ =
        (4.0/3.0)*sqr(phase_)*phase_.d()*g0_*(1.0 + e_)
       *sqrt(Theta_/constant::mathematical::pi);

       // 'thermal' conductivity
    kappa_ = conductivityModel_->kappa
    (
        phase_,
        Theta_,
        g0_,
        phase_.rho(),
        phase_.d(),
        e_
    );

    nuFric_ = frictionalStressModel_->nu
    (
        phase_,
        alphaMinFriction_,
        alphaMax_,
        frictionalStressModel_->frictionalPressure
        (
            phase_,
            alphaMinFriction_,
            alphaMax_
        )/phase_.rho(),
        symm(fvc::grad(phase_.U()))
    );

    // Limit viscosity and add frictional viscosity
    nu_.min(maxNut_);
    nuFric_ = min(nuFric_, maxNut_ - nu_);

    if (debug)
    {
        Info<< "    max(nu) = " << max(nu_).value() << nl
            << "    max(nuFric) = " << max(nuFric_).value() << endl;
    }
}

// ************************************************************************* //
