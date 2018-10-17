#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating twoPhaseSystem\n" << endl;

twoPhaseSystem fluid(mesh, g);

phaseModel& phase1 = fluid.phase1();
phaseModel& phase2 = fluid.phase2();

volScalarField& alpha1 = phase1;
volVectorField& U1 = phase1.U();
volVectorField& U2 = phase2.U();

volScalarField& p = phase1.thermo().p();

// Pressure
dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.U()
);

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
mesh.setFluxRequired(alpha1.name());

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
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

volScalarField rAU1
(
    IOobject
    (
        IOobject::groupName("rAU", phase1.name()),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimTime/dimDensity, 0.0)
);

volScalarField rAU2
(
    IOobject
    (
        IOobject::groupName("rAU", phase2.name()),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimTime/dimDensity, 0.0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K1(IOobject::groupName("K", phase1.name()), 0.5*magSqr(U1));
volScalarField K2(IOobject::groupName("K", phase2.name()), 0.5*magSqr(U2));

#include "createMRF.H"

fluid.pPrimeByA() = tmp<surfaceScalarField>
(
    new surfaceScalarField
    (
        IOobject
        (
            "pPrimeByA",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar
        (
            "0",
            dimensionSet(0, 2, -1, 0, 0),
            0.0
        )
    )
);
