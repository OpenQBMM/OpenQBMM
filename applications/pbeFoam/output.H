    runTime.write();

    volScalarField moment0(mesh.lookupObject<volScalarField>("moment.0"));

    Info<< "moment.0 = " << moment0[0] << endl;


