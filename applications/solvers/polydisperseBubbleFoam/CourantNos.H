#include "CourantNo.H"

{
    // Relative Courant number
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();
    scalar U1CoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();
    CoNum = max(max(CoNum, UrCoNum), U1CoNum);
}

{
    // relative advection Courant number
    scalar UiCoNum = 0.0;
    for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
    {
        UiCoNum =
            max
            (
                UiCoNum,
                0.5*gMax
                (
                    fvc::surfaceSum
                    (
                        mag(fvc::flux(phase1.Us(nodei)) - phi2)
                    )().primitiveField()/mesh.V().field()
                )*runTime.deltaTValue()
            );
    }
    CoNum = max(CoNum, UiCoNum);

    Info<< "Max U." << phase1.name() << " Courant Number = " << UiCoNum << endl;
}
