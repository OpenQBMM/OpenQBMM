    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

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

    Switch timeDependentVelocity
        = populationBalanceProperties.lookupOrDefault<bool>
            ("timeDependentVelocity", false);

    if (timeDependentVelocity)
    {
        #include "timeDependentVelocity.H"
    }

    Switch TaylorGreenVelocity
        = populationBalanceProperties.lookupOrDefault<bool>
            ("TaylorGreenVelocity", false);

    if (TaylorGreenVelocity)
    {
        #include "TaylorGreenVelocity.H"
    }

    autoPtr<populationBalanceModel> populationBalance
    (
        populationBalanceModel::New
        (
            "populationBalance", populationBalanceProperties, phi
        )
    );
