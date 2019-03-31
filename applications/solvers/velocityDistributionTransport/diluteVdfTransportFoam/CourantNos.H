scalar CoNum = 0.5*gMax
(
    fvc::surfaceSum
    (
        mag(fluid.phi()/fvc::interpolate(rho1))
    )().primitiveField()/mesh.V().field()
)*runTime.deltaTValue();

// Compute Courant number based on node velocities
CoNum = max(CoNum, dilutePhase->CoNum());

// Set max Courant number based on scheme
scalar pbeCoNum = dilutePhase->realizableCo();
maxCo = min(maxCo, pbeCoNum);

Info<< "Pbe Realizable Courant Number = " << pbeCoNum << endl;
Info<< "Max Courant Number = " << CoNum << endl;
