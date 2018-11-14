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

Switch solveP_rgh
(
    pimple.dict().lookupOrDefault<Switch>
    (
        "solveP_rgh", false
    )
);

label pRefCell = 0;
scalar pRefValue = 0.0;
if (solveP_rgh)
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );
    mesh.setFluxRequired(p_rgh.name());
}
else
{
    setRefCell
    (
        p,
        pimple.dict(),
        pRefCell,
        pRefValue
    );
    mesh.setFluxRequired(p.name());
}
mesh.setFluxRequired(alpha1.name());

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
