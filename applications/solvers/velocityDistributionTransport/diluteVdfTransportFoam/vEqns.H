for (label nodei = 0; nodei < phase2->nNodes(); nodei++)
{
    volScalarField alphaRho(phase2->volumeFraction(nodei)*rho2);
    const volVectorField& Ui = phase2->U(nodei);

    fvVectorMatrix UEqn
    (
        alphaRho*fvm::ddt(Ui)
      - alphaRho*fvc::ddt(Ui)
     ==
        AEqns[nodei]
    );

    UEqn.relax();
    UEqn.solve();
}

phase2->correct();
K2 = 0.5*magSqr(U2);
