// Solve the Momentum equation
MRF.correctBoundaryVelocity(U1);

fvVectorMatrix UEqn
(
    fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
  + MRF.DDt(alpha1*rho1 + Cvm, U1)
  + phase1->turbulence().divDevRhoReff(U1)
 ==
    virtualMassSource
  - Cvm
   *(
        fvm::ddt(U1)
      + fvm::div(phi1, U1)
      - fvm::Sp(fvc::div(phi1), U1)
    )
  + alpha1*rho1*g
  + dragSource
  - fvm::Sp(Kd, U1)
  + F
  + fvOptions(alpha1, rho1, U1)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
      - alpha1*fvc::grad(p)
    );

    fvOptions.correct(U1);
    K1 = 0.5*magSqr(U1);
}
