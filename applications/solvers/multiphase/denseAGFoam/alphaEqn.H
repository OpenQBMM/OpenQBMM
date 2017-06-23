{
    word alphaScheme("div(phi," + alpha1.name() + ')');
    const dictionary& alphaControls = mesh.solverDict
    (
        alpha1.name()
    );

    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    surfaceScalarField alpha1f("alpha1f", fvc::interpolate(alpha1));
/*
    surfaceScalarField alphaPhic1
    (
        fvc::flux
        (
            h2Fnf*phipS,
            alphap,
            alphaScheme
        )
    );

    explicitSolve(alphap,alphaPhic1,mesh.time().deltaTValue());
*/


    fvScalarMatrix alphaEqn(alpha1, alpha1.dimensions()*dimVol/dimTime);
    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {

		alpha1f = fvc::interpolate(max(alpha1, scalar(0)));

        alphaEqn =
        (
            fvm::ddt(alpha1)
          - fvc::ddt(alphap)
          + fvm::div(phi1, alpha1, alphaScheme)
          - fvm::laplacian(alpha1f*fluid.pPrimeByA()(), alpha1)
        );

        alphaEqn.relax();
        alphaEqn.solve();
    }

    phase1.alphaPhi() = alphaEqn.flux();

   	alpha1.max(0);

    alpha2 = scalar(1) - alpha1;
    phase2.alphaPhi() = fvc::interpolate(alpha2)*phi2;
    phase2.alphaRhoPhi() = phase2.alphaPhi()*fvc::interpolate(rho2);

    Info<< alpha1.name()
    << "  Average = "<< alpha1.weightedAverage(mesh.V()).value()
    << "  Min = " << min(alpha1).value()
    << "  Max = " << max(alpha1).value()
    << endl;

    ddtAlpha1Dilute = fvc::ddt(alpha1);

}
