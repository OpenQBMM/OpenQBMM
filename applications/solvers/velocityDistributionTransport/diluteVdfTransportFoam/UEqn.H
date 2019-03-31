// Solve the Momentum equation

MRF.correctBoundaryVelocity(U1);
fvVectorMatrix U1Eqn
(
    fvm::ddt(alpha1, rho1, U1)
  + fvm::div(alphaRhoPhi1, U1)
  + MRF.DDt(alpha1*rho1 + Cvm, U1)
  + fluid.turbulence().divDevRhoReff(U1)
 ==
    virtualMassSource
  - Cvm
   *(
        fvm::ddt(U1)
      + fvm::div(phi1, U1)
      - fvm::Sp(fvc::div(phi1), U1)
    )
  + dragSource
  - fvm::Sp(Kd, U1)
  + F
  + fvOptions(alpha1, rho1, U1)
);

U1Eqn.relax();
fvOptions.constrain(U1Eqn);
