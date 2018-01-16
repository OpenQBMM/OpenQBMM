scalar CoNum = 0;

// relative advection Courant number
forAll(U, nodei)
{
    CoNum =
        max
        (
            CoNum,
            0.5*gMax
            (
                fvc::surfaceSum
                (
                    mag
                    (
                        fvc::flux
                        (
                            U[nodei]*weights[nodei]
                           /(weights[nodei] + SMALL)
                        )
                    )
                )().primitiveField()/mesh.V().field()
            )*runTime.deltaTValue()
        );
}
Info<< "Max Courant Number = " << CoNum << endl;
