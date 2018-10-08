#include "CourantNo.H"

{
    // Relative Courant number
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    CoNum = max(CoNum, UrCoNum);

    scalar particleCoNumKin = (AGmodel.maxUxDx()*runTime.deltaTValue());
    scalar particleCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    Info << "Max particle phi Courant Number = " <<  particleCoNum  << endl;
    Info << "Max particle kinetic Courant Number = " <<  particleCoNumKin  << endl;
    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    CoNum = max(particleCoNum, max(particleCoNumKin, CoNum));
}
