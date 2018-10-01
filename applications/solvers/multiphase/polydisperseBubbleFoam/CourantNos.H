#include "CourantNo.H"

{
    // Relative Courant number
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    CoNum = max(CoNum, UrCoNum);
    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    // relative advection Courant number
    if (phase1.nNodes() > 1)
    {
        scalar maxNodeCoNum = 0;
        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            maxNodeCoNum =
                max
                (
                    maxNodeCoNum,
                    0.5*gMax
                    (
                        fvc::surfaceSum
                        (
                            mag(fvc::flux(phase1.Us(nodei)))
                        )().primitiveField()/mesh.V().field()
                    )*runTime.deltaTValue()
                );
        }
        CoNum = max(CoNum, maxNodeCoNum);
        Info<< "Max Node Courant Number = " << maxNodeCoNum << endl;
    }




}
