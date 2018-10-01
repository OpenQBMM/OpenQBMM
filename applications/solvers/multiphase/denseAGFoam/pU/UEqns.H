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

    volSymmTensorField Sigma2
    (
        IOobject::groupName("Sigma", phase2.name()),
        alpha2*nu2Eff*dev(twoSymm(fvc::grad(U2)))
    );

    U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1)
          + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(contErr1, U1)
          + MRF.DDt(alpha1*rho1, U1)

          + fluid.divDevRhoReff1()
         ==
            alpha1*rho2*fvc::div(Sigma2)
          + fvOptions(alpha1, rho1, U1)
        );

    U1Eqn.relax();
    U1Eqn += fvm::Sp(Kd, U1);
    fvOptions.constrain(U1Eqn);
    U1.correctBoundaryConditions();
    fvOptions.correct(U1);

    U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2)
          + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(contErr2, U2)
          + MRF.DDt(alpha2*rho2, U2)

          + fluid.divDevRhoReff2()
         ==
            fvOptions(alpha2, rho2, U2)
        );

    U2Eqn.relax();
    U2Eqn += fvm::Sp(Kd, U2);
    fvOptions.constrain(U2Eqn);
    U2.correctBoundaryConditions();
    fvOptions.correct(U2);
}
