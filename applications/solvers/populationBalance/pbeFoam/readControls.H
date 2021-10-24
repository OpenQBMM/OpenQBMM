    if (runTime.controlDict().lookupOrDefault("suppressSolverInfo", false))
    {
        lduMatrix::debug = 0;
    }
