forAll(nodes, nodei)
{
    const volVelocityNode& node = nodes[nodei];
    volScalarField alphaRho(node.primaryWeight()*rhop);
    volVectorField& Upi = nodes[nodei].velocityAbscissae();


    fvVectorMatrix UpiEqn
    (
        alphaRho*fvm::ddt(Upi)
      - alphaRho*fvc::ddt(Upi)
     ==
        g*alphaRho

      + Kds[nodei]*U
      - fvm::Sp(Kds[nodei], Upi)
    );

    UpiEqn.relax();
    UpiEqn.solve();
}

// Update moments to include drag and gravity
quadrature.updateMoments();
