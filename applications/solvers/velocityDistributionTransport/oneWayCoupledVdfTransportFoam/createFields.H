#include "readGravitationalAcceleration.H"

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

word dilutePhaseName(phaseProperties.lookup("dilutePhase"));
autoPtr<phaseModel> dilutePhase
(
    new vdfPhaseModel
    (
        mesh,
        phaseProperties,
        dilutePhaseName
    )
);

word fluidName(phaseProperties.lookup("continuousPhase"));
phaseModel fluid
(
    mesh,
    phaseProperties,
    fluidName
);
const volVectorField& U = fluid.U();

autoPtr<phasePair> pair
(
    new orderedPhasePair
    (
        dilutePhase,
        fluid,
        g,
        dimensionedScalar::lookupOrDefault
        (
            "sigma",
            phaseProperties,
            dimensionSet(1, 0, -2, 0, 0),
            0.0
        ),
        phaseProperties.subDict("aspectRatio")
    )
);
autoPtr<dragModel> drag
(
    dragModel::New
    (
        phaseProperties.subDict("drag"),
        pair()
    )
);
