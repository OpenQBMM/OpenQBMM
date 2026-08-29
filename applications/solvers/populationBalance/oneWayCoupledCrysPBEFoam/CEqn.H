{

    volScalarField k
    (
        IOobject("k_crys", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless/dimTime, 0.0)
    );

    volScalarField::Internal Su
    (
        IOobject("Su_crys", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimMass/dimVolume/dimTime, 0.0)
    );

    
    if (populationBalance->speciesCoupled())
    {
        const volScalarField::Internal& SYact = populationBalance->SYact();    
        const scalar smallC = 1.0e-12;
        forAll(k, celli)
        {
            scalar currentC = C[celli];
            scalar rate = SYact[celli];
            if (rate > 0.0 && currentC > smallC)
            {
                k[celli] = rate / currentC;
                Su[celli] = 0.0;
            }
            else if (rate < 0.0)
            {
                k[celli] = 0.0;
                Su[celli] = -rate; 
            }

            else
            {
                k[celli] = 0.0;
                Su[celli] = 0.0;
            }
        }
    }

    alphat = turbulence->nut()/Sct + mD;
    alphat.correctBoundaryConditions();

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi, C)
      - fvm::laplacian(alphat, C)
     ==
        Su
      - fvm::Sp(k, C)
      + fvOptions(C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);

    CEqn.solve();

    fvOptions.correct(C);

    C.correctBoundaryConditions();

    volScalarField Csat_field = saturationModel->Csat(T);
    dimensionedScalar minC("minC", dimDensity, 0.0);
    dimensionedScalar maxC("maxC", dimDensity, 20.0*gMax(Csat_field)); 
    
    C.clamp_range(minC, maxC);

    sigma = (C - Csat_field) / max(Csat_field, dimensionedScalar("minCsat", dimDensity, SMALL));
    sigma.correctBoundaryConditions();

    Info<< "C Equation solved. Max implicit rate constant k = " << gMax(k) << " s^-1. "
        << "Min/Max C = " << gMin(C) << " / " << gMax(C) << endl;
}
