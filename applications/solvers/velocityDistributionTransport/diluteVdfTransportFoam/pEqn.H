// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi1*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField alpharAUf("alpharAUf", fvc::interpolate(alpha1*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U1, p_rgh));

surfaceScalarField phig
(
    alpharAUf
   *(
       -ghf*fvc::snGrad(rho)*mesh.magSf()
      + fvc::interpolate(alpha2*(rho1 - rho2))*(g & mesh.Sf())
    )
);

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(alpha1*HbyA)
      + MRF.zeroFilter(alpharAUf*fvc::ddtCorr(alphaRho1, U1, alphaRhoPhi1))
    )
  + phig
);

MRF.makeRelative(fvc::interpolate(alphaRho1), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, alphaRho1, U1, phiHbyA, alpharAUf, MRF);

// fvScalarMatrix p_rghDDtEqn
// (
//     fvc::ddt(alphaRho1) + psi1*correction(fvm::ddt(p_rgh))
//   + fvc::div(phiHbyA)
//   ==
//     fvOptions(psi1, p_rgh, rho.name())
// );

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvc::div(phiHbyA)
      - fvm::laplacian(alpharAUf, p_rgh)
      + fvc::ddt(alpha1)
    );

    if (thermo1.incompressible())
    {
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
    }

    p_rghEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi1 = phiHbyA + p_rghEqn.flux();

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U1 = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/alpharAUf);
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
        K1 = 0.5*magSqr(U1);
    }
}

p = p_rgh + rho*gh;

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
      + fvc::div(alphaPhi, rho1)
      ==
        fvOptions(alpha1, rho)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho1);
}

alphaPhi1 = fvc::interpolate(alpha1)*phi1;
alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;

#include "contErrs.H"


