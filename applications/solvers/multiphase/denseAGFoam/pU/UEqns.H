Info<< "Constructing momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

volScalarField Kd(fluid.Kd());
{


    volScalarField nu2Eff(phase2.turbulence().nuEff());
//     if(cluterInducedViscosity)
//     {
//         nu2Eff += phase1.turbulence().nu()*(pow(alpha2, -2.8) - 1.0);
//     }

    nu2Eff *= alpha2;

    volTensorField gradU2(fvc::grad(U2));
    volSymmTensorField Sigma2
    (
        "Sigma2",
        nu2Eff*(twoSymm(gradU2) - (2.0/3.0)*tr(gradU2)*I)
    );
    volTensorField R2("R2",nu2Eff*gradU2 - Sigma2);

    U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(contErr1, U1)
          + MRF.DDt(alpha1*rho1, U1)
          + fluid.divDevRhoReff1()
         ==
            rho2*alpha1*fvc::div(Sigma2)
          + fvOptions(alpha1, rho1, U1)
        );

    U1Eqn.relax();
    U1Eqn += fvm::Sp(Kd, U1);
    fvOptions.constrain(U1Eqn);


    U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2)
          + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(contErr2, U2)
          + MRF.DDt(alpha2*rho2, U2)

          - fvm::laplacian(alpha2*rho2*nu2Eff, U2)
          + fvc::div(alpha2*rho2*R2)
         ==
            fvOptions(alpha2, rho2, U2)
        );

    U2Eqn.relax();
    U2Eqn += fvm::Sp(Kd, U2);
    fvOptions.constrain(U2Eqn);
}
