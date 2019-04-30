// Solve the Momentum equation
volScalarField alphaRho1(alpha1*rho1);
MRF.correctBoundaryVelocity(U1);

fvVectorMatrix UEqn
(
    fvm::ddt(alphaRho1, U1) + fvm::div(alphaRhoPhi1, U1)
  + MRF.DDt(alphaRho1 + Cvm, U1)
  + phase1->turbulence().divDevRhoReff(U1)
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

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
        ==
        fvc::reconstruct
        (
            (
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            )*fvc::interpolate(alpha1)*mesh.magSf()
          + fvc::interpolate(alpha2*(rho1 - rho2))*(g & mesh.Sf())
        )
    );

    fvOptions.correct(U1);
    K1 = 0.5*magSqr(U1);
}
