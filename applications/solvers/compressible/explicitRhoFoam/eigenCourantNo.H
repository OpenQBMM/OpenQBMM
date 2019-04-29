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
}

Info<< "Mean and max Courant Numbers = "
    << meanCoNum << ", " << CoNum << endl;
