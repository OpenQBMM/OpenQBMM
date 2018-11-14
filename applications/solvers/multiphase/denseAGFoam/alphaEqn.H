const dictionary& alphaControls = mesh.solverDict
(
    alpha1.name()
);

tmp<surfaceScalarField> pPrimeByA(h2f*fluid.pPrimeByA()());

label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
word alphaScheme("div(" + phi1.name() + ',' + alpha1.name() + ')');

alpha1.correctBoundaryConditions();
phase1.correctInflowOutflow(alphaPhi1);

fvScalarMatrix alpha1Eqn(alpha1, alpha1.dimensions()*dimVol/dimTime);

if (pPrimeByA.valid())
{
    surfaceScalarField phiP
    (
        pPrimeByA()*fvc::snGrad(alpha1, "bounded")*mesh.magSf()
    );
    surfaceScalarField phiS(h2f*phi1 - phiP);

    for (label acorr = 0; acorr < nAlphaCorr; acorr++)
    {
		surfaceScalarField alpha1f(fvc::interpolate(max(alpha1, scalar(0))));

        alpha1Eqn =
        (
            fvm::ddt(alpha1)
	      + fvm::div(phiS, alpha1, alphaScheme)
          - fvm::laplacian(alpha1f*pPrimeByA(), alpha1)
        );

        alpha1Eqn.relax();
        alpha1Eqn.solve();
    }
}
else
{
    alpha1Eqn =
    (
        fvm::ddt(alpha1)
      + fvm::div(phi1, alpha1)
    );
    alpha1Eqn.relax();
    alpha1Eqn.solve();

}

alphaPhi1 = alpha1Eqn.flux();
alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;

alphaPhi2 = phi - alphaPhi1;
phase2.correctInflowOutflow(alphaPhi2);
alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

Info<< alpha1.name() << " volume fraction = "
    << alpha1.weightedAverage(mesh.V()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;

// Ensure the phase-fractions are bounded
alpha1.max(0);
alpha1.min(1);

alpha2 = scalar(1) - alpha1;
alpha2.correctBoundaryConditions();

