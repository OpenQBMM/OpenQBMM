/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "BlendedInterfacialModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "bubblePressureModel.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvc.H"
#include "fvm.H"
#include "fixedValueFvsPatchFields.H"
#include "blendingMethod.H"
#include "HashPtrTable.H"
#include "UniformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::twoPhaseSystem
(
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    ),

    g_(g),

    phase2_
    (
        phaseModel::New
        (
            *this,
            *this,
            wordList(lookup("phases"))[1]
        )
    ),

    phase1_
    (
        phaseModel::New
        (
            *this,
            *this,
            wordList(lookup("phases"))[0]
        )
    ),

    nNodes_(phase1_->nNodes()),

    dgdt_
    (
        IOobject
        (
            "dgdt",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dgdt", dimless/dimTime, 0)
    )
{
    if (phase2_->nNodes() != 1)
    {
        FatalErrorInFunction
            << "Phase 2 is not monodisperse. Only one polydisperse phase" << nl
            << " can currently be handled and should be phase1."
            << exit(FatalError);
    }
    phi_ = calcPhi();

    phase2_().volScalarField::operator=(scalar(1) - phase1_());


    // Blending
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().dict().dictName(),
            blendingMethod::New
            (
                iter().dict(),
                wordList(lookup("phases"))
            )
        );
    }


    // Pairs

    phasePair::scalarTable sigmaTable(lookup("sigma"));
    phasePair::dictTable aspectRatioTable(lookup("aspectRatio"));

    pair_.set
    (
        new phasePair
        (
            phase1_(),
            phase2_(),
            g,
            sigmaTable
        )
    );

    pair1In2_.set
    (
        new orderedPhasePair
        (
            phase1_(),
            phase2_(),
            g,
            sigmaTable,
            aspectRatioTable
        )
    );

    pair2In1_.set
    (
        new orderedPhasePair
        (
            phase2_(),
            phase1_(),
            g,
            sigmaTable,
            aspectRatioTable
        )
    );


    // Models

    drag_.set
    (
        new BlendedInterfacialModel<dragModel>
        (
            lookup("drag"),
            (
                blendingMethods_.found("drag")
              ? blendingMethods_["drag"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_,
            false // Do not zero drag coefficent at fixed-flux BCs
        )
    );

    virtualMass_.set
    (
        new BlendedInterfacialModel<virtualMassModel>
        (
            lookup("virtualMass"),
            (
                blendingMethods_.found("virtualMass")
              ? blendingMethods_["virtualMass"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    lift_.set
    (
        new BlendedInterfacialModel<liftModel>
        (
            lookup("lift"),
            (
                blendingMethods_.found("lift")
              ? blendingMethods_["lift"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    wallLubrication_.set
    (
        new BlendedInterfacialModel<wallLubricationModel>
        (
            lookup("wallLubrication"),
            (
                blendingMethods_.found("wallLubrication")
              ? blendingMethods_["wallLubrication"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    turbulentDispersion_.set
    (
        new BlendedInterfacialModel<turbulentDispersionModel>
        (
            lookup("turbulentDispersion"),
            (
                blendingMethods_.found("turbulentDispersion")
              ? blendingMethods_["turbulentDispersion"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    bubblePressure_.set
    (
        new BlendedInterfacialModel<bubblePressureModel>
        (
            lookup("bubblePressure"),
            (
                blendingMethods_.found("bubblePressure")
              ? blendingMethods_["bubblePressure"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    phase1_->setModels();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::~twoPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::rho() const
{
    return
        phase1_()*phase1_->thermo().rho()
      + phase2_()*phase2_->thermo().rho();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::U() const
{
    return phase1_()*phase1_->U() + phase2_()*phase2_->U();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::calcPhi() const
{
    return
        fvc::interpolate(phase1_())*phase1_->phi()
      + fvc::interpolate(phase2_())*phase2_->phi();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Kd(const label nodei) const
{
    return drag_->K(nodei, 0);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kd() const
{
    tmp<volScalarField> tKd
    (
        new volScalarField
        (
            IOobject
            (
                "Kd",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "Kd",
                dimDensity*dimViscosity/sqr(dimLength),
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        tKd.ref() += Kd(nodei);
    }
    return tKd;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Kdf(const label nodei) const
{
    return drag_->Kf(nodei, 0);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Kdf() const
{
    tmp<surfaceScalarField> tKdf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Kd",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "Kd",
                dimDensity*dimViscosity/sqr(dimLength),
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        tKdf.ref() += Kdf(nodei);
    }
    return tKdf;
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Vm(const label nodei) const
{
    return virtualMass_->K(nodei, 0);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Vm() const
{
    tmp<volScalarField> tVm
    (
        new volScalarField
        (
            IOobject
            (
                "Vm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "Vm",
                dimDensity,
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        tVm.ref() += Vm(nodei);
    }
    return tVm;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Vmf(const label nodei) const
{
    return virtualMass_->Kf(nodei, 0);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Vmf() const
{
    tmp<surfaceScalarField> tVmf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Vmf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "Vmf",
                dimDensity,
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        tVmf.ref() += Vmf(nodei);
    }
    return tVmf;
}


Foam::tmp<Foam::volVectorField>
Foam::twoPhaseSystem::F(const label nodei) const
{
    volVectorField DDtU1
    (
        fvc::ddt(phase1_->U())
      + fvc::div(phase1_->phi(), phase1_->U())
      - fvc::div(phase1_->phi())*phase1_->U()
    );
    volVectorField DDtUi
    (
        fvc::ddt(phase1_->U())
      + fvc::div(phase1_->phi(), phase1_->Us(nodei))
      - fvc::div(phase1_->phi())*phase1_->Us(nodei)
    );

    return
        lift_->F<vector>(nodei, 0)
      + wallLubrication_->F<vector>(nodei, 0)
      - bubblePressure_->F<vector>(nodei, 0)

      // Force due to deviation from mean velocity
      - Kd(nodei)*phase1_->Vs(nodei)
      + Vm(nodei)
       *(
            DDtU1
          - DDtUi
        );
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::F() const
{
    tmp<volVectorField> tF
    (
        new volVectorField
        (
            IOobject
            (
                "F",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector
            (
                "F",
                dimensionSet(1, -2, -2, 0, 0),
                Zero
            )
        )
    );
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        tF.ref() += F(nodei);
    }
    return tF;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Ff(const label nodei) const
{
    volVectorField DDtU1
    (
        fvc::ddt(phase1_->U())
      + fvc::div(phase1_->phi(), phase1_->U())
      - fvc::div(phase1_->phi())*phase1_->U()
    );
    volVectorField DDtUi
    (
        fvc::ddt(phase1_->U())
      + fvc::div(phase1_->phi(), phase1_->Us(nodei))
      - fvc::div(phase1_->phi())*phase1_->Us(nodei)
    );
    return
        lift_->Ff(nodei, 0)
      + wallLubrication_->Ff(nodei, 0)
      - bubblePressure_->Ff(nodei, 0)

      // Force due to deviation from mean velocity
      + fvc::flux
        (
          - Kd(nodei)*phase1_->Vs(nodei)
          + Vm(nodei)
           *(
                DDtU1
              - DDtUi
            )
        );
}

Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Ff() const
{
    tmp<surfaceScalarField> tFf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Ff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "Ff",
                dimensionSet(1, 0, -2, 0, 0),
                Zero
            )
        )
    );
    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        tFf.ref() += Ff(nodei);
    }
    return tFf;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::D(const label nodei) const
{
    return turbulentDispersion_->D(nodei, 0);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::D() const
{
    return turbulentDispersion_->D();
}


Foam::tmp<Foam::fvVectorMatrix> Foam::twoPhaseSystem::divDevRhoReff1()
{
    if (phase1_().BGviscosity())
    {
        volScalarField rhoNuEff1
        (
            "rhoNuEff1",
            phase1_()
           *phase1_().d()
           *mag(phase1_().U() - phase2_().U())
           *sqrt(phase1_()*phase2_())
           *(phase1_().rho() + phase2_().rho()*virtualMass(phase1_()).Cvm())
        );

        return
            fvc::div(rhoNuEff1*dev2(T(fvc::grad(phase1_().U()))))
          - fvm::laplacian(rhoNuEff1, phase1_().U());
    }
    else
    {
        volVectorField& U1(phase1_().U());

        return phase1().turbulence().divDevRhoReff(U1);
    }
}


Foam::tmp<Foam::fvVectorMatrix> Foam::twoPhaseSystem::divDevRhoReff2()
{
    if (phase2_().BGviscosity())
    {
        volScalarField rhoNuEff2
        (
            "rhoNuEff2",
            phase2_().rho()
           *phase2_()
           *(
                phase2_().nu()
              + phase1_()/max(phase2_(), phase2_().residualAlpha())
               *(
                    phase1_().rho()/phase2_().rho()
                  + virtualMass(phase1_()).Cvm()
                )*phase1_().d()
               *mag(phase1_().U() - phase2_().U())
               *sqrt(phase1_()*phase2_())
               *pos0(phase2_() - 0.1)
            )
        );

        return
            fvc::div(rhoNuEff2*dev2(T(fvc::grad(phase2_().U()))))
          - fvm::laplacian(rhoNuEff2, phase2_().U());
    }
    else
    {
        volVectorField& U2(phase2_().U());

        return phase2().turbulence().divDevRhoReff(U2);
    }
}


void Foam::twoPhaseSystem::solve()
{
    const Time& runTime = mesh_.time();

    volScalarField& alpha1 = phase1_();
    volScalarField& alpha2 = phase2_();

    const surfaceScalarField& phi1 = phase1_->phi();
    const surfaceScalarField& phi2 = phase2_->phi();

    const dictionary& alphaControls = mesh_.solverDict
    (
        alpha1.name()
    );

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();

    surfaceScalarField phic("phic", phi_);
    surfaceScalarField phir("phir", phi1 - phi2);

    tmp<surfaceScalarField> alpha1alpha2f;

    if (pPrimeByA_.valid())
    {
        alpha1alpha2f =
            fvc::interpolate(max(alpha1, scalar(0)))
           *fvc::interpolate(max(alpha2, scalar(0)));

        surfaceScalarField phiP
        (
            pPrimeByA_()*fvc::snGrad(alpha1, "bounded")*mesh_.magSf()
        );

        phir += phiP;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", dgdt_.dimensions(), 0.0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        forAll(dgdt_, celli)
        {
            if (dgdt_[celli] > 0.0)
            {
                Sp[celli] -= dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
                Su[celli] += dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
            }
            else if (dgdt_[celli] < 0.0)
            {
                Sp[celli] += dgdt_[celli]/max(alpha1[celli], 1e-4);
            }
        }

        surfaceScalarField alphaPhic1
        (
            fvc::flux
            (
                phic,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        phase1_->correctInflowOutflow(alphaPhic1);

        if (nAlphaSubCycles > 1)
        {
            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhic10(alphaPhic1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi_,
                    alphaPhic10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    UniformField<scalar>(phase1_->alphaMax()),
                    zeroField()
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase1_->alphaPhi() = alphaPhic10;
                }
                else
                {
                    phase1_->alphaPhi() += alphaPhic10;
                }
            }

            phase1_->alphaPhi() /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhic1,
                Sp,
                Su,
                UniformField<scalar>(phase1_->alphaMax()),
                zeroField()
            );

            phase1_->alphaPhi() = alphaPhic1;
        }

        if (pPrimeByA_.valid())
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alpha1alpha2f()*pPrimeByA_(), alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            phase1_->alphaPhi() += alpha1Eqn.flux();
        }

        phase1_->alphaRhoPhi() =
            fvc::interpolate(phase1_->rho())*phase1_->alphaPhi();

        phase2_->alphaPhi() = phi_ - phase1_->alphaPhi();
        phase2_->correctInflowOutflow(phase2_->alphaPhi());
        phase2_->alphaRhoPhi() =
            fvc::interpolate(phase2_->rho())*phase2_->alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        // Ensure the phase-fractions are bounded
        alpha1.max(0);
        alpha1.min(1);

        alpha2 = scalar(1) - alpha1;
    }
}


void Foam::twoPhaseSystem::relativeTransport()
{
    if (nNodes_ > 1)
    {
        phase1_->relativeTransport();
    }
}


void Foam::twoPhaseSystem::averageTransport()
{
    PtrList<fvVectorMatrix> AEqns(nNodes_);

    if (nNodes_ == 1)
    {
        phase1_->averageTransport(AEqns);
        phase1_->correct();

        return;
    }

    // Liquid viscous stress
    volSymmTensorField taul(phase2_->turbulence().devRhoReff());

    // Acceleration of liquid phase
    volVectorField DDtU2
    (
        fvc::ddt(phase2_->U())
      + fvc::div(phase2_->phi(), phase2_->U())
      - fvc::div(phase2_->phi())*phase2_->U()
    );

    for (label nodei = 0; nodei < nNodes_; nodei++)
    {
        //  Build matrix to solve for velocity abscissae due to interfacial
        //  forces
        AEqns.set
        (
            nodei,
            new fvVectorMatrix
            (
                phase1_->Us(nodei),
                phase1_->Us(nodei).dimensions()*dimDensity*dimVol/dimTime
            )
        );

        const volScalarField& p(mesh_.lookupObject<volScalarField>("p"));
        volScalarField alphaRhoi(phase1_->alphas(nodei)*phase1_->rho());

        //  Implicit drag term added to velocity abscissae equations
        volScalarField Kd(this->Kd(nodei));

        // Interfacial forces
        AEqns[nodei] +=
            // Buoyancy
            g_*alphaRhoi
          + (
              - fvc::grad(p)
              + fvc::div(taul)
            )*phase1_->alphas(nodei)

            // Drag
          + Kd*phase2_->U()
          - fvm::Sp(Kd, phase1_->Us(nodei))

            // Virtual Mass
          + Vm(nodei)
           *(
                DDtU2
              - (
                    fvm::ddt(phase1_->Us(nodei))
                  + fvm::div(phase1_->phi(), phase1_->Us(nodei))
                  - fvm::Sp(fvc::div(phase1_->phi()), phase1_->Us(nodei))
                )
            )

            // Dispersion, lift, wall lubrication, and bubble pressure
          - turbulentDispersion_->F<vector>(nodei, 0)
          - lift_->F<vector>(nodei, 0)
          - wallLubrication_->F<vector>(nodei, 0)
          + bubblePressure_->F<vector>(nodei, 0);

    }

    phase1_->averageTransport(AEqns);
    phase1_->correct();

    phi_ = phase1_->alphaPhi() + phase2_->alphaPhi();
}


void Foam::twoPhaseSystem::correct()
{
    phase1_->correct();
    phase2_->correct();
}


void Foam::twoPhaseSystem::correctTurbulence()
{
    phase1_->turbulence().correct();
    phase2_->turbulence().correct();
}


bool Foam::twoPhaseSystem::read()
{
    bool readOK = regIOobject::read();

    bool readOK1 = phase1_->read(readOK);
    bool readOK2 = phase2_->read(readOK);

    return (readOK1 || readOK2);
}


const Foam::dragModel& Foam::twoPhaseSystem::drag(const phaseModel& phase) const
{
    return drag_->phaseModel(phase);
}


const Foam::virtualMassModel&
Foam::twoPhaseSystem::virtualMass(const phaseModel& phase) const
{
    return virtualMass_->phaseModel(phase);
}


const Foam::dimensionedVector& Foam::twoPhaseSystem::g() const
{
    return g_;
}


const Foam::dimensionedScalar& Foam::twoPhaseSystem::sigma() const
{
    return pair_->sigma();
}


// ************************************************************************* //
