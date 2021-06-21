PtrList<fvVectorMatrix> AEqns(phase2->nNodes());

volScalarField Kd
(
    IOobject
    (
        "Kd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimDensity/dimTime, 0.0)
);

volVectorField dragSource
(
    IOobject
    (
        "dragSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimAcceleration*dimDensity, Zero)
);

volScalarField Cvm
(
    IOobject
    (
        "Cvm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimDensity, 0.0)
);

volVectorField virtualMassSource
(
    IOobject
    (
        "virtualMassSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimAcceleration*dimDensity, Zero)
);

volVectorField F
(
    IOobject
    (
        "F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimAcceleration*dimDensity, Zero)
);

volScalarField Kh
(
    IOobject
    (
        "Kh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -3, -1, 0), 0.0)
);

{
    // Liquid viscous stress
    volSymmTensorField taul(phase1->turbulence().devRhoReff());

    // Acceleration of liquid phase
    volVectorField DDtU
    (
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    forAll(AEqns, nodei)
    {
        volScalarField alpha(phase2->volumeFraction(nodei));
        volScalarField alphaRho(alpha*rho2);

        const volVectorField& Ui = phase2->U(nodei);

        // WenYu drag
        AEqns.set
        (
            nodei,
            new fvVectorMatrix(Ui, Ui.dimensions()*dimDensity*dimVol/dimTime)
        );
        volScalarField Kdi(drag->K(nodei, 0));

        dragSource += Kdi*Ui;
        Kd += Kdi;

        // Interfacial forces
        AEqns[nodei] +=
            // Buoyancy
            (
               - fvc::grad(p)
               + rho2*g
               + fvc::div(taul)
            )*alpha

            // Drag
          + Kdi*U1 - fvm::Sp(Kdi, Ui);

        if (virtualMass.valid())
        {
            volScalarField Cvmi(virtualMass->K(nodei, 0));
            virtualMassSource +=
                Cvmi
               *(
                    fvc::ddt(Ui)
                  + fvc::div(fvc::flux(Ui), Ui)
                  - fvc::div(Ui)*Ui
                );
            Cvm += Cvmi;

            AEqns[nodei] +=
                Cvmi
               *(
                    DDtU
                  - (
                        fvm::ddt(Ui)
                      + fvm::div(fvc::flux(Ui), Ui)
                      - fvm::Sp(fvc::div(Ui), Ui)
                    )
                );
        }

        volVectorField Fi
        (
            IOobject
            (
                "Fi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimAcceleration*dimDensity, Zero)
        );

        if (turbulentDispersion.valid())
        {
            Fi += turbulentDispersion->F(nodei, 0);
        }
        if (lift.valid())
        {
            Fi += lift->F(nodei, 0);
        }
        if (wallLubrication.valid())
        {
            Fi += wallLubrication->F(nodei, 0);
        }

        F += Fi;
        AEqns[nodei] -= Fi;

        if (heatTransfer.valid())
        {
            Kh += heatTransfer->K(nodei, 0);
        }
    }
}
