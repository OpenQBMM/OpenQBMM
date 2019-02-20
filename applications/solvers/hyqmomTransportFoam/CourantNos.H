// Compute Courant number based on node velocities
scalar CoNum = populationBalance->CoNum();

// Set max Courant number based on scheme
scalar pbeCoNum = populationBalance->realizableCo();
maxCo = min(maxCo, pbeCoNum);

Info<< "Pbe Realizable Courant Number = " << pbeCoNum << endl;
Info<< "Max Courant Number = " << CoNum << endl;
