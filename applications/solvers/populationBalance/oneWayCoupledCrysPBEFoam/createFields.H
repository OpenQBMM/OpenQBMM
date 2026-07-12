#include "createRDeltaT.H"


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

Switch solveFlow
(
    pimple.dict().getOrDefault<bool>("solveFlow", true)
);



singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Expose nu as a registered volScalarField so aggregation sub-models that
// query mesh.lookupObject<volScalarField>("nu") work in reduced-order cases.
volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->nu()
);

#include "createMRF.H"
#include "createFvOptions.H"


Info<< "Reading field Temperature Concentration\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// kinematic turbulent thermal thermal conductivity m2/s
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Read populationBalanceProperties dictionary
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

// Get crystalCoeffs subdictionary
const dictionary& crystalCoeffs = 
    populationBalanceProperties.subDict
    (
        word(populationBalanceProperties.lookup("populationBalanceModel")) + "Coeffs"
    );

// Read transport properties from crystalCoeffs
dimensionedScalar Pr
(
    "Pr",
    dimless,
    crystalCoeffs.getOrDefault<scalar>("Pr", 1.0)
);

dimensionedScalar Prt
(
    "Prt",
    dimless,
    crystalCoeffs.getOrDefault<scalar>("Prt", 1.0)
);

dimensionedScalar mD
(
    "molarDiff",
    dimViscosity,
    crystalCoeffs.getOrDefault<scalar>("molarDiff", 1.0e-5)
);

dimensionedScalar Sct
(
    "Sct",
    dimless,
    crystalCoeffs.getOrDefault<scalar>("Sct", 0.7)
);


dimensionedScalar unitConverter(dimDensity, 10);
dimensionedScalar unitConverter2(inv(dimDensity),0.0094);

autoPtr<solutionSaturationModel> saturationModel
(
    solutionSaturationModel::New
    (
        crystalCoeffs.subDict("saturationModel"),
        mesh
    )
);

volScalarField Csat(saturationModel->Csat(T));

volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (C - Csat) / max(Csat, dimensionedScalar("minCsat", dimDensity, SMALL)),
    fvPatchFieldBase::zeroGradientType()
);
sigma.correctBoundaryConditions();

// Read crystal particle density from crystalCoeffs
uniformDimensionedScalarField rhop
(
    IOobject
    (
        "rhop",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dimensionedScalar
    (
        "rhop",
        dimDensity,
        crystalCoeffs.getOrDefault<scalar>("rhop", 1000.0)
    )
);

// Store rhop in the object registry so it can be found by populationBalance
rhop.store();

autoPtr<PDFTransportModels::populationBalanceModels::crystalPopulationBalance>
populationBalance
(
    new PDFTransportModels::populationBalanceModels::crystalPopulationBalance
    (
        "populationBalance",
        crystalCoeffs,
        phi
    )
);

