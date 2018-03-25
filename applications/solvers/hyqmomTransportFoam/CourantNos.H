// Compute Courant number based on node velocities
scalar CoNum = populationBalance->CoNum();

// Set max Courant number based on scheme
maxCo = min(maxCo, populationBalance->realizableCo());

Info<< "Max Courant Number = " << CoNum << endl;
