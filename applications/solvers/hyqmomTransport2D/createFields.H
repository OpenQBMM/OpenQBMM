// Pressure
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

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p.name());

Info<< "Reading quadratureProperties\n" << endl;
IOdictionary quadratureProperties
(
    IOobject
    (
        "quadratureProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
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

Switch collisions(populationBalanceProperties.lookup("collisions"));
scalar e(readScalar(populationBalanceProperties.lookup("e")));
dimensionedScalar tau(populationBalanceProperties.lookup("tau"));
dimensionedScalar rTau = 1.0/tau;

hyperbolicConditionalMomentInversion momentInverter
(
    quadratureProperties, 2
);

labelListList indicies(9, labelList(2,0));
indicies[0] = {1, 1};
indicies[1] = {1, 2};
indicies[2] = {1, 3};
indicies[3] = {2, 1};
indicies[4] = {2, 2};
indicies[5] = {2, 3};
indicies[6] = {3, 1};
indicies[7] = {3, 2};
indicies[8] = {3, 3};

labelListList mIndicies(10, labelList(2,0));
mIndicies[0] = {0, 0};
mIndicies[1] = {1, 0};
mIndicies[2] = {0, 1};
mIndicies[3] = {2, 0};
mIndicies[4] = {1, 1};
mIndicies[5] = {0, 2};
mIndicies[6] = {3, 0};
mIndicies[7] = {0, 3};
mIndicies[8] = {4, 0};
mIndicies[9] = {0, 4};

PtrList<volScalarField> weights(9);
PtrList<volVectorField> U(9);
PtrList<volScalarField> moments(10);

PtrList<surfaceScalarField> wOwn(9);
PtrList<surfaceScalarField> wNei(9);
PtrList<surfaceVectorField> UOwn(9);
PtrList<surfaceVectorField> UNei(9);

PtrList<volScalarField> fluxes(10);

forAll(indicies, nodei)
{
    label i = indicies[nodei][0];
    label j = indicies[nodei][1];

    word name = Foam::name(i) + Foam::name(j);
    weights.set
    (
        nodei,
        new volScalarField
        (
            IOobject
            (
                "weight." + name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("w" + name, dimless, 0.0)
        )
    );
    wOwn.set
    (
        nodei,
        new surfaceScalarField
        (
            IOobject
            (
                "weightOwn."+ name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("wOwn" + name, dimless, 0.0)
        )
    );
    wNei.set
    (
        nodei,
        new surfaceScalarField
        (
            IOobject
            (
                "weightNei."+ name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("wNei" + name, dimless, 0.0)
        )
    );

    U.set
    (
        nodei,
        new volVectorField
        (
            IOobject
            (
                "U." + name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("U." + name, dimVelocity, Zero)
        )
    );

    UOwn.set
    (
        nodei,
        new surfaceVectorField
        (
            IOobject
            (
                "UOwn."+ name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("UOwn." + name, dimVelocity, Zero)
        )
    );
    UNei.set
    (
        nodei,
        new surfaceVectorField
        (
            IOobject
            (
                "UNei."+ name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("UNei." + name, dimVelocity, Zero)
        )
    );
}

forAll(mIndicies, mi)
{
    label i = mIndicies[mi][0];
    label j = mIndicies[mi][1];
    word name = Foam::name(i) + Foam::name(j);

    moments.set
    (
        mi,
        new volScalarField
        (
            IOobject
            (
                "moment." + name,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh/*,
            dimensionedScalar
            (
                "moment" + name,
                pow(dimVelocity, i + j),
                0.0
            ),
            weights[0].boundaryField().types()*/
        )
    );

    fluxes.set
    (
        mi,
        new volScalarField
        (
            IOobject
            (
                "flux." + name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "flux" + name,
                pow(dimVelocity, i + j)/dimTime,
                0.0
            )
        )
    );
}
volScalarField Theta("Theta", (moments[3] + moments[5])/3.0);

#include "invertMoments.H"

surfaceScalarField own
(
    IOobject
    (
        "own",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("own", dimless, 1.0)
);

surfaceScalarField nei
(
    IOobject
    (
        "nei",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("nei", dimless, -1.0)
);
