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
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word dilutePhaseName(phaseProperties.lookup("dilutePhase"));
word fluidName(phaseProperties.lookup("continuousPhase"));
phaseModel fluid
(
    mesh,
    phaseProperties,
    fluidName
);
rhoThermo& thermo1 = fluid.thermo();
volScalarField& p = thermo1.p();
volScalarField& alpha1 = fluid;
volVectorField& U1 = fluid.U();
surfaceScalarField& phi1 = fluid.phi();
surfaceScalarField& alphaPhi1 = fluid.alphaPhi();
surfaceScalarField& alphaRhoPhi1 = fluid.alphaRhoPhi();
volScalarField& rho1 = thermo1.rho();
const volScalarField psi1 = thermo1.psi();

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K1
(
    IOobject::groupName("K", fluid.name()),
    0.5*magSqr(U1)
);

autoPtr<phaseModel> dilutePhase
(
    new vdfPhaseModel
    (
        mesh,
        phaseProperties,
        dilutePhaseName
    )
);
volScalarField& alpha2(dilutePhase());
rhoThermo& thermo2 = dilutePhase->thermo();
const volScalarField& rho2 = thermo2.rho();
const volVectorField& U2 = dilutePhase->U();
surfaceScalarField& phi2 = dilutePhase->phi();
surfaceScalarField& alphaPhi2 = dilutePhase->alphaPhi();
surfaceScalarField& alphaRhoPhi2 = dilutePhase->alphaRhoPhi();

volScalarField K2
(
    IOobject::groupName("K", dilutePhase->name()),
    0.5*magSqr(U2)
);

alpha1 = 1.0 - alpha2;
alphaPhi1 = fvc::interpolate(alpha1)*phi1;
alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;

alphaPhi2 = fvc::interpolate(alpha2)*phi2;
alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
p_rgh = p - rho*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;

if (thermo1.incompressible())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );
}

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    phaseProperties
);

mesh.setFluxRequired(p_rgh.name());

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
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

autoPtr<phasePair> pair
(
    new orderedPhasePair
    (
        dilutePhase,
        fluid,
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

#include "createMRF.H"
#include "createFvOptions.H"
