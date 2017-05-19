// Calculate momentum exchange coefficients

// Calculate drag coefficients for each bubble size
volScalarField Kd
(
    IOobject
    (
        "Kd",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
        "Kd",
        dimDensity*dimViscosity/sqr(dimLength),
        0.0
    )
);

//  Coefficient of virtual mass
volScalarField Cvm
(
    IOobject
    (
        "Cvm",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
        "Cvm",
        dimDensity,
        0
    )
);

//  Lift and wall lubrication force and correction term to account for
//  polydispersity
volVectorField F
(
    IOobject
    (
        "F",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector
    (
        "zero",
        dimForce/dimVolume,
        Zero
    )
);


//  Matrix for updating Us based on force relative to the mean bubble
//  velocity (drag and virtual mass)
PtrList<fvVectorMatrix> AEqns(nNodes);
