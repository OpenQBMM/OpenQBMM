const dictionary& alphaControls = mesh.solverDict
(
    alpha1.name()
);

volScalarField& dgdt = fluid.dgdt();
surfaceScalarField pPrimeByA(h2f*fluid.pPrimeByA());

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

word alphaScheme("div(phi," + alpha1.name() + ')');
word alpharScheme("div(phir," + alpha1.name() + ')');

alpha1.correctBoundaryConditions();

surfaceScalarField phic("phic", phi);
surfaceScalarField phir("phir", phi1 - phi2);

tmp<surfaceScalarField> alpha1alpha2f;

if (fluid.pPrimeByA().valid())
{
    alpha1alpha2f =
        fvc::interpolate(max(alpha1, scalar(0)))
       *fvc::interpolate(max(alpha2, scalar(0)));

    surfaceScalarField phiP
    (
        fluid.pPrimeByA()()*fvc::snGrad(alpha1, "bounded")*mesh.magSf()
    );

    phir += phiP;
}

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
        dimensionedScalar("Sp", dgdt.dimensions(), 0.0)
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

    forAll(dgdt, celli)
    {
        if (dgdt[celli] > 0.0)
        {
            Sp[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
            Su[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        }
        else if (dgdt[celli] < 0.0)
        {
            Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
        }
    }

    // Flux in standard form with the addition of subtraction of dilute flux
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
          - fvc::flux(-phir, scalar(1) - alpha1, alpharScheme)
          - phi1*(1.0 - h2f),
            alpha1,
            alpharScheme
        )
    );

    phase1.correctInflowOutflow(alphaPhic1);

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
                phase1.alphaMax(),
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
            phase1.alphaMax(),
            0
        );

        alphaPhi1 = alphaPhic1;
    }

    if (fluid.pPrimeByA().valid())
    {
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1) - fvc::ddt(alpha1)
          - AGmodel.ddtAlphaDilute()
          - fvm::laplacian
            (
                alpha1alpha2f()*pPrimeByA,
                alpha1,
                "bounded"
            )
        );

        alpha1Eqn.relax();
        alpha1Eqn.solve();

        alphaPhi1 += alpha1Eqn.flux();
    }

    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;

    alphaPhi2 = phi - alphaPhi1;
    phase2.correctInflowOutflow(alphaPhi2);
    alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

    Info<< alpha1.name() << " volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    // Ensure the phase-fractions are bounded
    alpha1.max(0);
    alpha1.min(1);

    alpha2 = scalar(1) - alpha1;
}
