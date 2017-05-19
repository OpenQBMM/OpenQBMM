    Info<< "Calculating field DDtU1 and DDtU2\n" << endl;

    fvVectorMatrix DDtU1
    (
        fvm::ddt(U1)
      + fvm::div(phi1, U1)
      - fvm::Sp(fvc::div(phi1),U1)
    );

    fvVectorMatrix DDtU2
    (
        fvm::ddt(U2)
      + fvm::div(phi2, U2)
      - fvm::Sp(fvc::div(phi2),U2)
    );
