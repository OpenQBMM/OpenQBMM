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
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

volVectorField Up
(
    IOobject
    (
        "U.particles",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, Zero)
);

volSymmTensorField Sigmap
(
    IOobject
    (
        "Sigma.particles",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("Sigma", sqr(dimVelocity), Zero)
);

volScalarField Thetap
(
    IOobject
    (
        "Theta.particles",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Theta", sqr(dimVelocity), Zero)
);

autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "particles", populationBalanceProperties, phi
    )
);

velocityQuadratureApproximation& quadrature =
    mesh.lookupObjectRef<velocityQuadratureApproximation>
    (
        "quadratureProperties.particles"
    );

const labelListList& momentOrders = quadrature.momentOrders();
volVelocityMomentFieldSet& moments = quadrature.moments();
mappedPtrList<volVelocityNode>& nodes = quadrature.nodes();
const labelList& velocityIndexes(nodes[0].velocityIndexes());
label sizeIndex = nodes[0].sizeIndex();

bool computeVariance = false;

forAll(momentOrders, mi)
{
    forAll(velocityIndexes, cmpt)
    {
        if (momentOrders[mi][velocityIndexes[cmpt]] >= 2)
        {
            computeVariance = true;
            Sigmap.writeOpt() = IOobject::AUTO_WRITE;
            Thetap.writeOpt() = IOobject::AUTO_WRITE;
        }
    }
}

tmp<volScalarField> dMean;

if (sizeIndex != -1)
{
    dMean = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "d.particles",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("d", dimLength, 0.0)
        )
    );
}

#include "computeParticleFields.H"
