    if (mesh.nCells() != 1)
    {
        FatalErrorIn(args.executable())
            << "Solver only applicable to single cell cases"
            << exit(FatalError);
    }

    Info<< nl << "Reading thermophysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero),
        p.boundaryField().types()
    );

    #include "createPhi.H"

    Info << "Creating turbulence model.\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Reading populationBalanceProperties\n" << endl;

    IOdictionary populationBalanceProperties
    (
        IOobject
        (
            "populationBalanceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    autoPtr<populationBalanceModel> populationBalance
    (
        populationBalanceModel::New
        (
            "populationBalance",
            populationBalanceProperties,
            phi
        )
    );
