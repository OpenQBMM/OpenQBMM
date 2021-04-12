#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Reading phaseProperties\n" << endl;

IOdictionary phaseProperties
(
    IOobject
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

word dilutePhaseName(phaseProperties.lookup("dilutePhase"));
word fluidName(phaseProperties.lookup("continuousPhase"));

autoPtr<phaseModel> phase1
(
    new phaseModel
    (
        mesh,
        phaseProperties,
        fluidName
    )
);

rhoThermo& thermo1 = phase1->thermo();
volScalarField& p = thermo1.p();
volScalarField& alpha1 = phase1();
volVectorField& U1 = phase1->U();
surfaceScalarField& phi1 = phase1->phi();
surfaceScalarField& alphaPhi1 = phase1->alphaPhi();
surfaceScalarField& alphaRhoPhi1 = phase1->alphaRhoPhi();
volScalarField& rho1 = thermo1.rho();
const volScalarField psi1 = thermo1.psi();

Info<< "Creating field kinetic energy K\n" << endl;

volScalarField K1
(
    IOobject::groupName("K", phase1->name()),
    0.5*magSqr(U1)
);

autoPtr<phaseModel> phase2
(
    new vdfPhaseModel
    (
        mesh,
        phaseProperties,
        dilutePhaseName
    )
);

volScalarField& alpha2(phase2());
rhoThermo& thermo2 = phase2->thermo();
volScalarField& rho2 = thermo2.rho();
const volScalarField psi2 = thermo2.psi();
volVectorField& U2 = phase2->U();
surfaceScalarField& phi2 = phase2->phi();
surfaceScalarField& alphaPhi2 = phase2->alphaPhi();
surfaceScalarField& alphaRhoPhi2 = phase2->alphaRhoPhi();

volScalarField K2
(
    IOobject::groupName("K", phase2->name()),
    0.5*magSqr(U2)
);

alpha1 = 1.0 - alpha2;
alphaPhi1 = fvc::interpolate(alpha1)*phi1;
alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;

alphaPhi2 = fvc::interpolate(alpha2)*phi2;
alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);

mesh.setFluxRequired(p.name());

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    phaseProperties
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("zeroDpDt", p.dimensions()/dimTime, 0.0)
);

autoPtr<phasePair> pair
(
    new orderedPhasePair
    (
        phase2,
        phase1,
        g,
        dimensionedScalar::lookupOrDefault
        (
            "sigma",
            phaseProperties,
            dimensionSet(1, 0, -2, 0, 0),
            0.0
        ),
        phaseProperties.subDict("aspectRatio")
    )
);

autoPtr<dragModel> drag
(
    dragModel::New
    (
        phaseProperties.subDict("drag"),
        pair()
    )
);

autoPtr<virtualMassModel> virtualMass;

if (phaseProperties.found("virtualMass"))
{
    virtualMass =
        virtualMassModel::New
        (
            phaseProperties.subDict("virtualMass"),
            pair()
        );
}

autoPtr<liftModel> lift;

if (phaseProperties.found("lift"))
{
    lift =
        liftModel::New
        (
            phaseProperties.subDict("lift"),
            pair()
        );
}

autoPtr<wallLubricationModel> wallLubrication;

if (phaseProperties.found("wallLubrication"))
{
    wallLubrication =
        wallLubricationModel::New
        (
            phaseProperties.subDict("wallLubrication"),
            pair()
        );
}

autoPtr<turbulentDispersionModel> turbulentDispersion;

if (phaseProperties.found("turbulentDispersion"))
{
    turbulentDispersion =
        turbulentDispersionModel::New
        (
            phaseProperties.subDict("turbulentDispersion"),
            pair()
        );
}

autoPtr<heatTransferModel> heatTransfer;

if (phaseProperties.found("heatTransfer"))
{
    heatTransfer =
        heatTransferModel::New
        (
            phaseProperties.subDict("heatTransfer"),
            pair()
        );
}

#include "createMRF.H"
#include "createFvOptions.H"
