// Compute Courant number based on node velocities
CoNum = 0.0;
CoNum = max(CoNum, dilutePhase->CoNum());

// Set max Courant number based on scheme
scalar pbeCoNum = dilutePhase->realizableCo();
maxCo = min(maxCo, pbeCoNum);

Info<< "Pbe Realizable Courant Number = " << pbeCoNum << endl;
Info<< "Max Courant Number = " << CoNum << endl;
