scalar CoNum = 0.5*gMax
(
    fvc::surfaceSum
    (
        mag(phi1/fvc::interpolate(rho1))
    )().primitiveField()/mesh.V().field()
)*runTime.deltaTValue();

// Compute Courant number based on node velocities
CoNum = max(CoNum, phase2->CoNum());

// Set max Courant number based on scheme
scalar pbeCoNum = phase2->realizableCo();
maxCo = min(maxCo, pbeCoNum);

Info<< "Pbe Realizable Courant Number = " << pbeCoNum << endl;
Info<< "Max Courant Number = " << CoNum << endl;
