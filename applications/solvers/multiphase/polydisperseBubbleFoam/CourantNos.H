#include "CourantNo.H"

{
    // Relative Courant number
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    // relative advection Courant number
    for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
    {
        UrCoNum =
            max
            (
                UrCoNum,
                0.5*gMax
                (
                    fvc::surfaceSum
                    (
                        mag(fvc::flux(phase1.Us(nodei)) - phi2)
                    )().primitiveField()/mesh.V().field()
                )*runTime.deltaTValue()
            );
    }

    CoNum = max(CoNum, UrCoNum);

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;
}
