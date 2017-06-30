Info<< "Constructing momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

volScalarField Kd(fluid.Kd());
// volScalarField Kdi("Kdi", fluid.drag(phase1).Ki(0,0));
// volScalarField CdRe("CdRe", fluid.drag(phase1).CdRe(0,0));
// if (runTime.outputTime())
// {
//     Kd.write();
//     Kdi.write();
//     CdRe.write();
// }
{


    volScalarField nu2Eff(phase2.turbulence().nuEff());
//     if(cluterInducedViscosity)
//     {
//         nu2Eff += phase1.turbulence().nu()*(pow(alpha2, -2.8) - 1.0);
//     }

    volSymmTensorField Sigma2
    (
        "Sigma2",
        nu2Eff*dev(twoSymm(fvc::grad(U2)))
    );

    U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1)
          + fvm::div
            (
                AGmodel.hydrodynamicScalef(alphaRhoPhi1),
                U1,
                "div(" + alphaRhoPhi1.name() + "," + U1.name() +")"
            )
          - fvm::Sp(contErr1, U1)
          + MRF.DDt(alpha1*rho1, U1)

          + fluid.divDevRhoReff1()
         ==
            alpha1*rho2/rho1*fvc::div(phase1.turbulence().devRhoReff())
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

          + fluid.divDevRhoReff2()
         ==
            fvOptions(alpha2, rho2, U2)
        );

    U2Eqn.relax();
    U2Eqn += fvm::Sp(Kd, U2);
    fvOptions.constrain(U2Eqn);
}
