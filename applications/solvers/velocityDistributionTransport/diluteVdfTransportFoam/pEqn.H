// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi1*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField alpharAUf("alpharAUf", fvc::interpolate(alpha1*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U1, p));

surfaceScalarField phiHbyA1
(
    "phiHbyA1",
    fvc::flux(HbyA)
);

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(alpha1)*phiHbyA1
);

MRF.makeRelative(fvc::interpolate(alphaRho1), phiHbyA);

surfaceScalarField rAUf
(
    "rAUf",
    mag
    (
        fvc::interpolate(alpha1)*alpharAUf
    )
);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, alphaRho1, U1, phiHbyA, alpharAUf, MRF);

fvScalarMatrix pEqnComp1
(
    (
        contErr1
      - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
    )/rho1
  + (alpha1*psi1/rho1)*correction(fvm::ddt(p))
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvc::div(phiHbyA)
      - fvm::laplacian(rAUf, p)
      + fvc::ddt(alpha1)
      + pEqnComp1
    );

    if (thermo1.incompressible())
    {
        pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));
    }

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        surfaceScalarField mSfGradp("mSfGradp", pEqn.flux()/rAUf);

        // Calculate the conservative fluxes
        phi1 = phiHbyA1 + alpharAUf*mSfGradp;

        // Explicitly relax pressure for momentum corrector
        p.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U1 =
            HbyA
          + fvc::reconstruct
            (
                alpharAUf*(mSfGradp)
            );
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
        K1 = 0.5*magSqr(U1);
    }
}

// Thermodynamic density update
thermo1.correctRho(psi1*p - psip0);

if (thermo1.dpdt())
{
    dpdt = fvc::ddt(p);
}

if (!thermo1.incompressible())
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(alpha1, rho1)
      + fvc::div(alphaPhi1, rho1)
      ==
        fvOptions(alpha1, rho1)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho1);
}

alphaPhi1 = fvc::interpolate(alpha1)*phi1;
alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;

#include "contErrs.H"


