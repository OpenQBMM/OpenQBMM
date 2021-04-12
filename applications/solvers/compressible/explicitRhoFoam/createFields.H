#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;

compressibleSystem fluid(mesh);

fluidThermo& thermo = fluid.thermo();

const volScalarField& rho = fluid.rho();

const volVectorField& U = fluid.U();
const surfaceScalarField& phi = fluid.phi();

autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        fluid.massFlux(),
        fluid.thermo()
    )
);

const volScalarField& muEff(fluid.thermo().mu());
bool inviscid = (max(muEff).value() < SMALL);

errorEstimator error(mesh);
