#include "CourantNo.H"

{
    // Relative Courant number
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();
    CoNum = max(CoNum, UrCoNum);

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;
}

{
    //- mean advection Courant number
    scalar U1CoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    for (label nodei = 0; nodei < nNodes; nodei++)
    {
        U1CoNum =
            max
            (
                U1CoNum,
                0.5*gMax
                (
                    fvc::surfaceSum
                    (
                        mag(fvc::flux(phase1.Us(nodei)))
                    )().primitiveField()/mesh.V().field()
                )*runTime.deltaTValue()
            );
    }
    CoNum = max(CoNum, U1CoNum);

    Info<< "Max U." << phase1.name() << " Courant Number = " << U1CoNum << endl;
}

{
    // relative advection Courant number
    scalar VCoNum = 0.0;
    for (label nodei = 0; nodei < nNodes; nodei++)
    {
        VCoNum =
            max
            (
                VCoNum,
                0.5*gMax
                (
                    fvc::surfaceSum
                    (
                        mag(fvc::flux(phase1.Vs(nodei)))
                    )().primitiveField()/mesh.V().field()
                )*runTime.deltaTValue()
            );
    }
    CoNum = max(CoNum, VCoNum);

    Info<< "Max V Courant Number = " << VCoNum << endl;
}
