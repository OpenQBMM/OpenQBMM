 {
    const dictionary& alphaControls = mesh.solverDict(alpha1.name());

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phi1 - phi2);

    tmp<surfaceScalarField> alpha1alpha2f;

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", inv(dimTime), 0.0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi)*min(alpha1, scalar(1))
        );

        surfaceScalarField alphaPhic1
        (
            fvc::flux
            (
                phic,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        phase1.correctInflowFlux(alphaPhic1);

        if (nAlphaSubCycles > 1)
        {
            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhic10(alphaPhic1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi,
                    alphaPhic10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    1.0,
                    0
                );

                if (alphaSubCycle.index() == 1)
                {
                    alphaPhi1 = alphaPhic10;
                }
                else
                {
                    alphaPhi1 += alphaPhic10;
                }
            }

           alphaPhi1 /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                alphaPhic1,
                Sp,
                Su,
                1.0,
                0
            );

            alphaPhi1 = alphaPhic1;
        }

        alphaRhoPhi1 =
            fvc::interpolate(rho1)*alphaPhi1;

        alphaPhi2 = phi - alphaPhi1;
        phase2.correctInflowFlux(alphaPhi2);
        alphaRhoPhi2 =
            fvc::interpolate(rho2)*alphaPhi2;

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        alpha1.max(0);
        alpha1.min(1);

        alpha2 = scalar(1) - alpha1;
    }
}
