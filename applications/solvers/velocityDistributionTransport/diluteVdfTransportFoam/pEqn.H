surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

volScalarField rAU1
(
    IOobject::groupName("rAU", phase1->name()),
    1.0/(U1Eqn.A())
);

surfaceScalarField alpharAUf1
(
    fvc::interpolate(alpha1*rAU1)
);

// --- Pressure corrector loop
while (pimple.correct())
{
    // Update continuity errors due to temperature changes
    #include "correctContErrs.H"

    rho = alpha1*rho1 + alpha2*rho2;

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF.correctBoundaryFlux(U1, phi1);

    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1->name()),
        U1
    );
    HbyA1 = rAU1*U1Eqn.H();

    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    surfaceScalarField phig1
    (
        alpharAUf1
       *(
           ghSnGradRho
         - alphaf2*fvc::interpolate(rho1 - rho2)*(g & mesh.Sf())
        )
    );


    // ddtPhiCorr filter -- only apply in pure(ish) phases
    surfaceScalarField alphaf1Bar(fvc::interpolate(fvc::average(alphaf1)));
    surfaceScalarField phiCorrCoeff1(pos0(alphaf1Bar - 0.99));

    {
        surfaceScalarField::Boundary& phiCorrCoeff1Bf =
            phiCorrCoeff1.boundaryFieldRef();

        forAll(mesh.boundary(), patchi)
        {
            // Set ddtPhiCorr to 0 on non-coupled boundaries
            if
            (
               !mesh.boundary()[patchi].coupled()
             || isA<cyclicAMIFvPatch>(mesh.boundary()[patchi])
            )
            {
                phiCorrCoeff1Bf[patchi] = 0;
            }
        }
    }

    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1->name()),
        fvc::flux(HbyA1)
      + phiCorrCoeff1*fvc::interpolate(alpha1.oldTime()*rho1.oldTime()*rAU1)
       *(
            MRF.absolute(phi1.oldTime())
          - fvc::flux(U1.oldTime())
        )/runTime.deltaT()
      + fvc::flux(rAU1*F)
      - phig1
    );

    // Construct the mean predicted flux
    // including explicit drag contributions based on absolute fluxes
    surfaceScalarField phiHbyA("phiHbyA", alphaf1*phiHbyA1);
    MRF.makeRelative(phiHbyA);

    // Construct pressure "diffusivity"
    surfaceScalarField rAUf("rAUf", mag(alphaf1*alpharAUf1));

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf1.boundaryField()*phi1.boundaryField()
              + alphaf2.boundaryField()*phi2.boundaryField()
            )
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    tmp<fvScalarMatrix> pEqnComp1;

    // Construct the compressibility parts of the pressure equation
    if (pimple.transonic())
    {
        surfaceScalarField phid1
        (
            IOobject::groupName("phid", phase1->name()),
            fvc::interpolate(psi1)*phi1
        );

        pEqnComp1 =
            (
                contErr1
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + correction
            (
                (alpha1/rho1)*
                (
                    psi1*fvm::ddt(p_rgh)
                  + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                )
            );
        pEqnComp1.ref().relax();
    }
    else
    {
        pEqnComp1 =
            (
                contErr1
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));
    }

    // Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple.correctNonOrthogonal())
    {
        // Construct the transport part of the pressure equation
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          + fvc::ddt(alpha1)
        );

        pEqnIncomp.setReference(pRefCell, pRefValue);

        solve
        (
            pEqnIncomp
        );

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi1 = phiHbyA1 + alphaf1*pEqnIncomp.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U1 = fvc::reconstruct(phi1);
            U1.correctBoundaryConditions();
            fvOptions.correct(U1);
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh + rho*gh, pMin);

    // Limit p_rgh
    p_rgh = p - rho*gh;

    // Update densities from change in p_rgh
    rho1 += psi1*(p_rgh - p_rgh_0);

    // Correct p_rgh for consistency with p and the updated densities
    rho = alpha1*rho1 + alpha2*rho2;
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();


}

// Update the phase kinetic energies
K1 = 0.5*magSqr(U1);

// Update the pressure time-derivative if required
if (thermo1.dpdt())
{
    dpdt = fvc::ddt(p);
}

alphaPhi1 = fvc::interpolate(alpha1)*phi1;
alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
