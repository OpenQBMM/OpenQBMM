for (label nodei = 0; nodei < dilutePhase->nNodes(); nodei++)
{
    const volVectorField& Ui = dilutePhase->U(nodei);
    volScalarField Kd(drag->Ki(nodei, 0)/dilutePhase->rho());

    fvVectorMatrix UEqn
    (
        fvm::ddt(Ui)
      - fvc::ddt(Ui)
     ==
        Kd*U
      - fvm::Sp(Kd, Ui)
    );

    UEqn.relax();
    UEqn.solve();
}
dilutePhase->correct();
