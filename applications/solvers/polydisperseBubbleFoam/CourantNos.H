#include "CourantNo.H"

{
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    for (label nodei = 0; nodei < nNodes; nodei++)
    {
        scalarField sumPhii
        (
            fvc::surfaceSum
            (
                mag(fvc::flux(phase1.Us(nodei)) - phi2)
            )().primitiveField()
        );

        UrCoNum =
            max
            (
                UrCoNum,
                0.5*gMax(sumPhii/mesh.V().field())*runTime.deltaTValue()
            );
    }
    CoNum = max(CoNum, UrCoNum);

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;
}
