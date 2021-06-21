// Solve the change in velocity momentums due to drag
for (label nodei = 0; nodei < dilutePhase->nNodes(); nodei++)
{
    volScalarField alphaRho
    (
      dilutePhase->volumeFraction(nodei)*dilutePhase->rho()
    );

    const volVectorField& Ui = dilutePhase->U(nodei);
    volScalarField Kd(drag->K(nodei, 0));

    fvVectorMatrix UEqn
    (
        alphaRho*fvm::ddt(Ui)
      - alphaRho*fvc::ddt(Ui)
     ==
        Kd*U
      - fvm::Sp(Kd, Ui)
      + alphaRho*g
    );

    UEqn.relax();
    UEqn.solve();
}

dilutePhase->correct();
