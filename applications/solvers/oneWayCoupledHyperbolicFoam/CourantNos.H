// Compute Courant number based on node velocities
CoNum = populationBalance->CoNum();

// Compute Courant number based on granular temperature
CoNum =
    max
    (
        CoNum,
        0.5*gMax
        (
            fvc::surfaceSum
            (
                mag(fvc::interpolate(sqrt(Thetap))*mesh.Sf())
            )().primitiveField()/mesh.V().field()
        )*mesh.time().deltaTValue()
    );

// Set max Courant number based on scheme
scalar pbeCoNum = populationBalance->realizableCo();
maxCo = min(maxCo, pbeCoNum);

Info<< "Pbe Realizable Courant Number = " << pbeCoNum << endl;
Info<< "Max Courant Number = " << CoNum << endl;
