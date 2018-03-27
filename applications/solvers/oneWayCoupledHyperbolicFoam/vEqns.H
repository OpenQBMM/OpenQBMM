forAll(nodes, nodei)
{
    volVectorField& Upi = nodes[nodei].primaryAbscissa();
    fvVectorMatrix UpiEqn
    (
        fvm::ddt(Upi)
      - fvc::ddt(Upi)
     ==
        1.0/St*U
      - fvm::Sp(1.0/St, Upi)
    );

    UpiEqn.relax();
    UpiEqn.solve();
}

//- Update moments to include drag source
quadrature.updateMoments();