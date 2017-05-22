phase1.correct();

// Set drag coefficient, Virtual mass coefficient, and F to Zero so that
//  individual size contributions can be summed
Kd ==
    dimensionedScalar
    (
        "0",
        Kd.dimensions(),
        0.0
    );

Cvm ==
    dimensionedScalar
    (
        "0",
        Cvm.dimensions(),
        0
    );

F ==
    dimensionedVector
    (
        "0",
        F.dimensions(),
        Zero
    );

//  Liquid viscous stress tensor (same for all sizes)
volTensorField taul
(
    "taul",
    rho2*alpha2*nu2*(fvc::grad(U2) + dev2(T(fvc::grad(U2))))
);

//  Set source terms
for (label nodei = 0; nodei < nNodes; nodei++)
{
    volScalarField alphaRhoi =
        max(phase1.alphas(nodei), phase1.residualAlpha())*rho1;

    AEqns.set
    (
        nodei,
        new fvVectorMatrix
        (
            phase1.Us(nodei),
            phase1.Us(nodei).dimensions()*dimVol/dimTime
        )
    );


    //  Buoyancy force for velocity moment transport
    {
        AEqns[nodei] +=
        (
            g
          + (
              - fvc::grad(p)
              + fvc::div(taul)
            )/rho1
        );
    }


    //  Drag force
    {
        volScalarField Kdi = drag->K(nodei,0);
        Kd += Kdi;

        //  Relative drag term, added as a source term to the two fluid
        F -= Kdi*phase1.Vs(nodei);

        //  Implicit drag term added to velocity abscissae equations
        AEqns[nodei] +=
        (
            Kdi/alphaRhoi*U2
          - fvm::Sp
            (
                Kdi/alphaRhoi,
                phase1.Us(nodei)
            )
        );
    }

    //  Disperson force
    //  Dispersion for mean velocity is handeld in pEqn.H
    {
        AEqns[nodei] +=
            dispersion->F<vector>(nodei,0)/alphaRhoi;
    }


    // Virtual Mass
    {
        volScalarField Cvmi = virtualMass->K(nodei,0);
        Cvm += Cvmi;

        surfaceScalarField phiAlpha = fvc::flux(phase1.Us(nodei));

        fvVectorMatrix DDtUs
        (
            fvm::ddt(phase1.Us(nodei))
          + fvm::div(phi1, phase1.Us(nodei))
          - fvm::Sp(fvc::div(phi1), phase1.Us(nodei))
        );

        volScalarField DDtMi
        (
            Cvmi
           *(
                fvc::ddt(phase1.alphas(nodei))
              + fvc::div(phi1, phase1.alphas(nodei))
              - fvc::div(phi1)*phase1.alphas(nodei)
            )
        );

        F +=
        (
            Cvmi
           *(
                (DDtU1 & U1)
              - (DDtUs & phase1.Us(nodei))
            )
          + phase1.Vs(nodei)*DDtMi
        );

        AEqns[nodei] +=
        (
            Cvmi/alphaRhoi
           *(
                (DDtU2 & U2)
              - DDtUs
            )
          + fvm::Sp(DDtMi/alphaRhoi, phase1.Us(nodei))
          - U2*DDtMi/alphaRhoi
        );
    }


    // Lift Force
    {
        volVectorField Fli = lift->F<vector>(nodei,0);

        F += Fli;
        AEqns[nodei] += Fli/alphaRhoi;
    }


    // Wall lubrication
    {
        volVectorField Fwli = wallLubrication->F<vector>(nodei,0);

        F += Fwli;
        AEqns[nodei] += Fwli/alphaRhoi;
    }
}
