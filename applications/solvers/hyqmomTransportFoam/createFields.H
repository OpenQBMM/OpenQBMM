//- Empty phi field to pass into poplation balance model
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("zero", dimVolume/dimTime, 0.0)
);

Info<< "Reading populationBalanceProperties\n" << endl;
IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
volScalarField dp
(
    IOobject
    (
        "d.particles",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("d", dimLength, 1e-4)
);
autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "particles", populationBalanceProperties, phi
    )
);


