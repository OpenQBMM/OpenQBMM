PtrList<fvVectorMatrix> AEqns(dilutePhase->nNodes());
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

{
    // Liquid viscous stress
    volSymmTensorField taul(fluid.turbulence().devRhoReff());

    // Acceleration of liquid phase
    volVectorField DDtU
    (
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    forAll(AEqns, nodei)
    {
        volScalarField alpha(dilutePhase->volumeFraction(nodei));
        volScalarField alphaRho(alpha*rho2);

        const volVectorField& Ui = dilutePhase->U(nodei);

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
               - gh*fvc::grad(rho) - alpha1*(rho2 - rho1)*g
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
            Fi -= turbulentDispersion->F(nodei, 0);
        }
        if (lift.valid())
        {
            Fi -= lift->F(nodei, 0);
        }
        if (wallLubrication.valid())
        {
            Fi -= wallLubrication->F(nodei, 0);
        }

        F += Fi;
        AEqns[nodei] -= Fi;
    }
}
