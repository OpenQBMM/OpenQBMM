{
    surfaceScalarField amaxSf
    (
        mag(phi) + fvc::interpolate(fluid.speedOfSound())*mesh.magSf()
    );

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    CoNum = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumAmaxSf)/gSum(mesh.V().field()))*runTime.deltaTValue();

    scalar realizableCo = populationBalance().realizableCo();

    Info<< "Max Courant Number for realizable moment transport  = "
        << realizableCo << endl;

     maxCo = min(maxCo, realizableCo);
}

Info<< "Mean and max Courant Numbers = "
    << meanCoNum << ", " << CoNum << endl;
