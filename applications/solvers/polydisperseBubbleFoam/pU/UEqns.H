Info<< "Constructing momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

{
    // Set U1 and U2 equations
    U1Eqn =
    (
        fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
      - fvm::Sp(contErr1, U1)

     ==

      - Cvm
       *(
            DDtU1
          - (DDtU2 & U2)
        )
      + fvc::div((alpha1*rho1*nu1)*dev2(T(fvc::grad(U1))))
      + fvm::laplacian(alpha1*rho1*nu1, U1)
    );
    U1Eqn.relax();
    U1Eqn += fvm::Sp(Kd, U1);
    U1.correctBoundaryConditions();


//     tmp<volScalarField> nuEff =
//         nu2 + alpha1/max(alpha2, phase2.residualAlpha())*(rho1/rho2)*nu1;

    U2Eqn =
    (
        fvm::ddt(alpha2, rho2, U2) + fvm::div(alphaRhoPhi2, U2)
      - fvm::Sp(contErr2, U2)

      ==

        fvc::div((alpha2*rho2*nu2)*dev2(T(fvc::grad(U2))))
      + fvm::laplacian(alpha2*rho2*nu2, U2)

      - Cvm
        *(
            DDtU2
          - (DDtU1 & U1)
        )
    );
    U2Eqn.relax();
    U2Eqn += fvm::Sp(Kd, U2);
    U2.correctBoundaryConditions();
}
