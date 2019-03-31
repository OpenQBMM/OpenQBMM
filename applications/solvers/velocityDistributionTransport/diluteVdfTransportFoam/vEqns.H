for (label nodei = 0; nodei < dilutePhase->nNodes(); nodei++)
{
    volScalarField alphaRho(dilutePhase->volumeFraction(nodei)*rho2);
    const volVectorField& Ui = dilutePhase->U(nodei);

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
dilutePhase->correct();
K2 = 0.5*magSqr(U2);
