{
    scalar time = runTime.value();

    volVectorField cellCentres(mesh.C());

    forAll(cellCentres, celli)
    {
        scalar x = cellCentres[celli].x();

        U[celli] = vector((1.0 - x)/(1.0 + time), 0, 0);
    }

    U.correctBoundaryConditions();
    U.write();
    phi = fvc::flux(U);
}
