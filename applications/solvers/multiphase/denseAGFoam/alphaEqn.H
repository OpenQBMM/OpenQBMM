
{
    word alphaScheme("div(phi," + alpha1.name() + ')');
    int nAlphaCorr(readInt(pimple.dict().lookup("nAlphaCorr")));

    surfaceScalarField alpha1f("alpha1f", fvc::interpolate(alpha1));

    tmp<surfaceScalarField> pPrimeByA =
        h2Fnf*fvc::interpolate(rAU1*phase1.turbulence().pPrime());
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

	phi1S *=  h2Fnf;

    fvScalarMatrix alphaEqn(alpha1, alpha1.dimensions()*dimVol/dimTime);
    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {

		alpha1f = fvc::interpolate(max(alpha1, scalar(0)));

        alphaEqn =
        (
            fvm::ddt(alpha1)
//			- fvc::ddt(alphap)
		    + fvm::div(phi1S, alpha1,alphaScheme)
            - fvm::laplacian(alpha1f*pPrimeByA(), alpha1)
        );

        alphaEqn.relax();
        alphaEqn.solve();
    }

//    alphaPhip = alphaPhic1 + alphapEqn.flux();
    phase1.alphaPhi() = alphaEqn.flux();

   	alpha1.max(0);
	alpha1f = fvc::interpolate(alpha1);
    phase1.alphaPhi() = alpha1f*phi1;
    phase1.alphaRhoPhi() = phase1.alphaPhi()*fvc::interpolate(rho1);

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
