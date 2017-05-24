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
        tmp<volScalarField> Kdi = drag->K(nodei,0);
        tmp<volScalarField> Ki = drag->Ki(nodei,0)/rho1;
        Kd += Kdi();

        //  Relative drag term, added as a source term to the two fluid
        F -= Kdi()*phase1.Vs(nodei);

        //  Implicit drag term added to velocity abscissae equations
        AEqns[nodei] +=
        (
            Ki()*U2
          - fvm::Sp
            (
                Ki(),
                phase1.Us(nodei)
            )
        );
    }

    //  Disperson force
    //  Dispersion for mean velocity is handeld in pEqn.H
    {
        AEqns[nodei] += dispersion->A(nodei,0);
    }


    // Virtual Mass
    {
        tmp<volScalarField> Cvmi = virtualMass->K(nodei,0);
        tmp<volScalarField> Ci = virtualMass->Ki(nodei,0)/rho1;

        Cvm += Cvmi();

        fvVectorMatrix DDtUs
        (
            fvm::ddt(phase1.Us(nodei))
          + fvm::div(phi1, phase1.Us(nodei))
          - fvm::Sp(fvc::div(phi1), phase1.Us(nodei))
        );

        volScalarField DDtMi
        (
            fvc::ddt(phase1.alphas(nodei))
          + fvc::div(phi1, phase1.alphas(nodei))
          - fvc::div(phi1)*phase1.alphas(nodei)
        );

        F +=
        (
            Cvmi()
           *(
                (DDtU1 & U1)
              - (DDtUs & phase1.Us(nodei))
              + phase1.Vs(nodei)*DDtMi
            )
        );

        volScalarField alphaRhoi
            (
                max
                (
                    phase1.alphas(nodei),
                    phase1.residualAlpha()
                )*rho1
            );

        AEqns[nodei] +=
        (
            Ci()
           *(
                (DDtU2 & U2)
              - DDtUs
            )
          + fvm::Sp
            (
                Ci()*DDtMi,
                phase1.Us(nodei)
            )
          - U2*Ci()*DDtMi
        );
    }


    // Lift Force
    {
        F += lift->F<vector>(nodei,0);
        AEqns[nodei] += lift->A(nodei,0);
    }


    // Wall lubrication
    {
        F += lift->F<vector>(nodei,0);
        AEqns[nodei] += lift->A(nodei,0);
    }
}
