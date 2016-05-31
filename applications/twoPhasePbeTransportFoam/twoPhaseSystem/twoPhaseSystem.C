/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
#include "heatTransferModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "fvMatrix.H"
#include "surfaceInterpolate.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcCurl.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fixedValueFvsPatchFields.H"

#include "blendingMethod.H"
#include "HashPtrTable.H"


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

    phase1_
    (
        *this,
        *this,
        wordList(lookup("phases"))[0]
    ),

    phase2_
    (
        *this,
        *this,
        wordList(lookup("phases"))[1]
    ),

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
        this->calcPhi()
    ),

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
    ),
    populationBalance_
    (
        Foam::populationBalanceModel::New
        (
            "populationBalance",
            mesh.lookupObject<IOdictionary>("populationBalanceProperties"),
            phase2_.U(),
            phi_
        )
    )
{
    phase2_.volScalarField::operator=(scalar(1) - phase1_);


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
            phase1_,
            phase2_,
            g,
            sigmaTable
        )
    );

    pair1In2_.set
    (
        new orderedPhasePair
        (
            phase1_,
            phase2_,
            g,
            sigmaTable,
            aspectRatioTable
        )
    );

    pair2In1_.set
    (
        new orderedPhasePair
        (
            phase2_,
            phase1_,
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

    heatTransfer_.set
    (
        new BlendedInterfacialModel<heatTransferModel>
        (
            lookup("heatTransfer"),
            (
                blendingMethods_.found("heatTransfer")
              ? blendingMethods_["heatTransfer"]
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::~twoPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::rho() const
{
    return phase1_*phase1_.thermo().rho() + phase2_*phase2_.thermo().rho();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::U() const
{
    return phase1_*phase1_.U() + phase2_*phase2_.U();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::calcPhi() const
{
    return
        fvc::interpolate(phase1_)*phase1_.phi()
      + fvc::interpolate(phase2_)*phase2_.phi();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kd() const
{
    return drag_->K();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Kdf() const
{
    return drag_->Kf();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Vm() const
{
    return virtualMass_->K();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Vmf() const
{
    return virtualMass_->Kf();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kh() const
{
    return heatTransfer_->K();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::F() const
{
    return lift_->F<vector>() + wallLubrication_->F<vector>();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Ff() const
{
    return lift_->Ff() + wallLubrication_->Ff();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::D() const
{
    return turbulentDispersion_->D();
}


void Foam::twoPhaseSystem::solve()
{
    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;
    dimensionedScalar one("one", alpha1.dimensions(), 1.0);
    
    populationBalance_->solve();
    
    const volScalarField& alphaGas =
        mesh_.lookupObject<volScalarField>("moment.1.populationBalance");
        
    alpha1 = alphaGas;
    alpha2 = one - alpha1;

    Info<< alpha1.name() << " volume fraction = "
        << alpha1.weightedAverage(mesh_.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}


void Foam::twoPhaseSystem::correct()
{
    phase1_.correct();
    phase2_.correct();
}


void Foam::twoPhaseSystem::correctTurbulence()
{
    phase1_.turbulence().correct();
    phase2_.turbulence().correct();
}


bool Foam::twoPhaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        readOK &= phase1_.read(*this);
        readOK &= phase2_.read(*this);

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
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


const Foam::dimensionedScalar& Foam::twoPhaseSystem::sigma() const
{
    return pair_->sigma();
}


// ************************************************************************* //
