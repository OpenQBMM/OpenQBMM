{
    volVectorField cellCentres(mesh.C());

    forAll(cellCentres, celli)
    {
        scalar x = cellCentres[celli].x();
        scalar y = cellCentres[celli].y();
        scalar twoPi = 2.0*Foam::constant::mathematical::pi;

        U[celli] =
            vector
            (
                Foam::sin(twoPi*x)*Foam::cos(twoPi*y),
                - Foam::cos(twoPi*x)*Foam::sin(twoPi*y),
                0
            );
    }

    U.correctBoundaryConditions();
    U.write();
    phi = fvc::flux(U);
}
