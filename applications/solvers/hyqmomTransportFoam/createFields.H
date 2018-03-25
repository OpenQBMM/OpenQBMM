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

Info<< "Reading phaseProperties\n" << endl;

// Default field used for Boltzmann collision
volScalarField dp
(
    IOobject
    (
        IOobject::groupName("d", "particles"),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("d", dimLength, 3e-4)
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

autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "particles", populationBalanceProperties, phi
    )
);


