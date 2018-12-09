while (pimple.correct())
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    rAU1 = 1.0/U1Eqn.A();
    rAU2 = 1.0/U2Eqn.A();

    surfaceScalarField alpharAUf1(fvc::interpolate(alpha1*rAU1));
    surfaceScalarField alpharAUf2(fvc::interpolate(alpha2*rAU2));

    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1.name()),
        rAU1*U1Eqn.H()
    );


    volVectorField HbyA2
    (
        IOobject::groupName("HbyA", phase2.name()),
        rAU2*U2Eqn.H()
    );

    // particle pressure flux (e.g. due to particle-particle pressure)
    if (implicitPhasePressure)
    {
        fluid.pPrimeByA() = fvc::interpolate(rAU1*phase1.turbulence().pPrime());
    }

    surfaceScalarField phiPp
    (
        "phiPp",
        fluid.pPrimeByA()()*fvc::snGrad(alpha1)*mesh.magSf()
    );

    forAll(U1.boundaryField(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>(U1.boundaryField()[patchi])
        || isA<wallFvPatch>(mesh.boundary()[patchi])
        )
        {
            phiPp.boundaryFieldRef()[patchi] == 0.0;
        }
    }

    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        (fvc::interpolate(HbyA1) & mesh.Sf())
    + alpharAUf1*fvc::interpolate(rho1)*fvc::ddtCorr(U1, phi1)
    );

    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        (fvc::interpolate(HbyA2) & mesh.Sf())
    + alpharAUf2*fvc::interpolate(rho2)*fvc::ddtCorr(U2, phi2)
    );

    volScalarField rAUKd1(rAU1*Kd);
    phiHbyA1 +=
    (
        fvc::interpolate(rAUKd1)*phi2
    - phiPp
    + alpharAUf1*fvc::interpolate(rho1)*(g & mesh.Sf())
    );

    volScalarField rAUKd2(rAU2*Kd);
    phiHbyA2 +=
    (
        fvc::interpolate(rAUKd2)*phi1
    + alpharAUf2*fvc::interpolate(rho2)*(g & mesh.Sf())
    );

    surfaceScalarField phiHbyA("phiHbyA", h2f*alpha1f*phiHbyA1 + alpha2f*phiHbyA2);

    HbyA1 += rAUKd1*U2;
    HbyA2 += rAUKd2*U1;

    surfaceScalarField rAUf
    (
        "rAUf",
        mag(h2f*alpha1f*alpharAUf1 + alpha2f*alpharAUf2)
    );


    // Update the pressure BCs to ensure flux consistency
    U = AGmodel.h2()*alpha1*U1 + alpha2*U2;
    constrainPressure(p, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p)
          + AGmodel.ddtAlphaDilute()
        );
        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp("mSfGradp", pEqn.flux()/rAUf);

            phi1.boundaryFieldRef() =
                    mesh.Sf().boundaryField() & U1.boundaryField();

            phi1 = phiHbyA1 + alpharAUf1*mSfGradp;

            phi2.boundaryFieldRef() =
                    mesh.Sf().boundaryField() & U2.boundaryField();

            phi2 = phiHbyA2 + alpharAUf2*mSfGradp;

            phi = alpha1f*phi1 + alpha2f*phi2 ;

            p.relax();
            mSfGradp = pEqn.flux()/rAUf;

            U1 = HbyA1
            + fvc::reconstruct
                (
                    alpharAUf1
                *(
                        (g & mesh.Sf())*fvc::interpolate(rho1)
                    + mSfGradp
                    )
                - phiPp
                );
            U1.correctBoundaryConditions();
            fvOptions.correct(U1);

            U2 = HbyA2
            + fvc::reconstruct
                (
                    alpharAUf2
                *(
                        (g & mesh.Sf())*fvc::interpolate(rho2)
                    + mSfGradp
                    )
                );
            U2.correctBoundaryConditions();
            fvOptions.correct(U2);

        }
    }
}
