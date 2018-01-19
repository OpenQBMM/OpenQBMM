// Mean velocities
volScalarField m00 = max(moments[0], SMALL);
volScalarField u("u", moments[1]/m00);
volScalarField uSqr("uSqr", sqr(u));
volScalarField v("v", moments[2]/m00);
volScalarField vSqr("vSqr", sqr(v));

// Variances of velocities
dimensionedScalar zeroVar("zeroVar", sqr(dimVelocity), 0.0);
volScalarField sigma1
(
    "sigma1",
    max(moments[3]/m00 - uSqr, zeroVar)
);
volScalarField sigma2
(
    "sigma2",
    max(moments[5]/m00 - vSqr, zeroVar)
);
Theta = (sigma1 + sigma2)/3.0;

scalar omega = (1.0 + e)/2.0;
scalar omegaSqr = sqr(omega);
scalar a1 = omegaSqr;
scalar b1 = a1 - 2.0*omega + 1.0;

volScalarField sigma11
(
    "sigma11",
    a1*Theta + b1*sigma1
);
volScalarField sigma22
(
    "sigma22",
    a1*Theta + b1*sigma2
);
volScalarField sigma12
(
    "sigma12",
    b1*(moments[4]/m00 - u*v)
);

PtrList<volScalarField> Meq(10);
forAll(mIndicies, mi)
{
    label i = mIndicies[mi][0];
    label j = mIndicies[mi][1];
    word name = Foam::name(i) + Foam::name(j);

    Meq.set
    (
        mi,
        new volScalarField
        (
            IOobject
            (
                "Meq." + name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "Meq" + name,
                pow(dimVelocity, i + j),
                0.0
            )
        )
    );
}

Meq[0] = moments[0];
Meq[1] = moments[1];
Meq[2] = moments[2];
Meq[3] = m00*(sigma11 + uSqr);
Meq[4] = m00*(sigma12 + u*v);
Meq[5] = m00*(sigma22 + vSqr);
Meq[6] = m00*(3.0*sigma11*u + u*uSqr);
Meq[7] = m00*(3.0*sigma22*v + v*vSqr);
Meq[8] = m00*(6.0*uSqr*sigma11 + 3.0*sqr(sigma11) + uSqr*uSqr);
Meq[8] = m00*(6.0*vSqr*sigma22 + 3.0*sqr(sigma22) + vSqr*vSqr);
