#include "readGravitationalAcceleration.H"

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

// Particle fields
Info<< "Reading phaseProperties\n" << endl;
IOdictionary phaseProperties
(
    IOobject
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
word particleName = "particles";//(phaseProperties.lookup("dilutePhase"));
scalar residualAlpha(phaseProperties.lookupType<scalar>("residualAlpha"));
dimensionedScalar St
(
    "St",
    dimTime,
    phaseProperties
);
Info<<particleName << endl;
Info<<phi.name()<<endl;
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
        particleName, populationBalanceProperties, phi
    )
);
velocityQuadratureApproximation& quadrature =
    mesh.lookupObjectRef<velocityQuadratureApproximation>
    (
        IOobject::groupName
        (
            "quadratureProperties",
            particleName
        )
    );

const labelListList& momentOrders = quadrature.momentOrders();
volVectorMomentFieldSet& moments = quadrature.moments();
mappedPtrList<volVectorNode>& nodes = quadrature.nodes();
const label nDims = momentOrders[0].size();

const volScalarField& alpha = moments(labelList(nDims, 0));

volVectorField Up
(
    IOobject
    (
        IOobject::groupName("U", particleName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, Zero)
);

volScalarField Thetap
(
    IOobject
    (
        IOobject::groupName("Theta", particleName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", sqr(dimVelocity), 0.0)
);

#include "computeParticleFields.H"
