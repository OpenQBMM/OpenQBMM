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

word particleName = "particles";
volVectorField Up
(
    IOobject
    (
        IOobject::groupName("U", particleName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
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
const volVelocityMomentFieldSet& moments = quadrature.moments();
mappedPtrList<volVelocityNode>& nodes = quadrature.nodes();
const labelList& velocityIndexes = nodes[0].velocityIndexes();
const label sizeIndex = nodes[0].sizeIndex();
const label nDims = momentOrders[0].size();

#include "computeParticleFields.H"
