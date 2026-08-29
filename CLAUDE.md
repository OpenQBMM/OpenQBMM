# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

OpenQBMM is an OpenFOAM **module** (not a standalone project) implementing Quadrature-Based
Moment Methods (QBMM) for polydisperse multiphase flows: population balance, turbulent mixing,
and velocity-distribution (kinetic) transport. It targets the **openfoam.com (ESI) line**;
`master` and `development-openfoam.com` are the live branches. CI (`.github/workflows/main.yml`)
builds against OpenFOAM 2412.

Everything requires a sourced OpenFOAM environment (`WM_PROJECT_DIR` must be set) — every build
script starts by sourcing `$WM_PROJECT_DIR/wmake/scripts/AllwmakeParseArguments`.

## Build

```sh
./Allwmake -j                 # libraries + applications
./Allwclean                   # wclean all of src, applications, test
```

Or, without a sourced environment (as CI does):

```sh
/usr/bin/openfoam2412 ./Allwmake -j
```

Artifacts go to `FOAM_MODULE_LIBBIN`/`FOAM_MODULE_APPBIN`, which default to
`FOAM_USER_LIBBIN`/`FOAM_USER_APPBIN`. Override with `-prefix` or the `FOAM_MODULE_*` env vars;
`FOAM_MODULE_PREFIX=false` skips the build entirely. Every `Make/options` repeats the
`sinclude $(GENERAL_RULES)/module-path-user` + failsafe block — copy that boilerplate into any
new `Make/options`.

`Allwmake` has `buildTests=false` hardcoded (disabled for releases), so tests must be built
explicitly. Rebuilding a single library or app: `wmake libso src/quadratureMethods/momentInversion`,
`wmake applications/solvers/populationBalance/pbeFoam`.

**Build order matters.** `src/Allwmake` lists the libraries in dependency order
(Vandermonde → momentSets → momentInversion → hermiteQuadrature → quadratureNode →
fieldMomentInversion → quadratureApproximations → momentAdvection → PDFTransportModels →
mixingModels → populationBalanceModels). A new library must be inserted at the right point,
and inter-library includes are written as relative paths into the other library's `lnInclude`
(e.g. `-I../momentSets/lnInclude`), so `lnInclude` for a dependency must exist before dependents
compile.

## Tests

```sh
./test/Allwmake -j     # build unit-test executables into FOAM_USER_APPBIN
./test/Allrun          # run the unit-test suite
```

Tests are OpenFOAM applications named `Test-*` (`Test-Vandermonde`, `Test-UnivariateMomentSet`,
`Test-UnivariateMomentInversion`, `Test-PopulationBalanceModel`, ...). Run one directly by name
once built. Tests that read dictionaries must be run from their own directory — e.g.

```sh
cd test/momentInversions/univariate/univariateMomentInversion && Test-UnivariateMomentInversion
```

which reads `quadraturePropertiesGauss`, `...GQMOM`, `...Lobatto`, `...Radau` from that directory.
Tests with a `testCase/` subdirectory (`populationBalanceModel`, `univariateQuadratureApproximation`)
are run from inside `testCase`. `test/univariateMomentAdvection` holds many mesh-refinement cases,
each with its own `Allrun`/`Allclean` — these are convergence studies, not part of `test/Allrun`.

`tutorials/` and `validation/` are per-solver OpenFOAM cases, each with `Allrun`/`Allclean`;
`validation/` cases additionally produce comparison graphs against analytical/reference data.
Note that `test/Allwmake` deliberately skips `conditionalMomentInversion` (commented out, needs
updating for the current structure).

## Architecture

The core abstraction chain, roughly bottom-up:

- **`Vandermonde`** — linear algebra for moment problems (Vandermonde system solution).
- **`momentSets`** (`univariateMomentSet`, `multivariateMomentSet`) — a set of moment *values*
  plus realizability checking. `supportType` (a scoped enum in `src/quadratureMethods/supportType`)
  records the support of the distribution (R, R+, [0,1], ...) and drives which inversion/realizability
  rules apply.
- **`momentInversion`** — the algorithms that turn moments into quadrature weights/abscissae.
  - `univariate/basic`: Gauss, Gauss-Radau, Gauss-Lobatto, generalized (GQMOM), hyperbolic (HyQMOM).
  - `univariate/extended`: EQMOM kernels — beta, gamma, lognormal.
  - `multivariate`: CHyQMOM, CHyQMOM+, sizeCHyQMOM, TensorProduct, monoKinetic, conditional.
  These operate on plain `momentSet`s — no mesh, no fields. This is what the unit tests exercise.
- **`fieldMomentInversion`** — applies a `momentInversion` cell-by-cell over an `fvMesh`
  (`basicFieldMomentInversion`, `basicVelocityFieldMomentInversion`). This is the bridge from
  pointwise algorithm to CFD field.
- **`quadratureNode`** / **`moments`** — the field-level containers: a node holds weight + abscissae
  (+ sigma for EQMOM) as `volScalarField`s / `surfaceScalarField`s; `momentFieldSet` holds the moment
  fields.
- **`quadratureApproximation`** (templated on `<momentType, nodeType>`) — the central object. It is
  an `IOdictionary` read from `constant/quadratureProperties[.<name>]` and owns the moment fields,
  the node list, and the `fieldMomentInversion`. Solvers and models talk to this, not to the
  inversion algorithms directly.
- **`momentAdvection`** — spatial transport of moments in a realizable way (`univariate`:
  firstOrderKinetic, zeta, none; `velocity`: kinetic-based schemes with their own `fvQuadraturePatch`
  boundary handling).
- **`PDFTransportModels`** — abstract base tying advection + source terms together
  (`univariatePDFTransportModel`, `velocityPDFTransportModel`).
- **`populationBalanceModels`** / **`mixingModels`** — the concrete physics on top of PDF transport,
  selected at run time (`univariate`, `velocity`, `sizeVelocity`, `mixingPopulationBalance`, `none`;
  and `turbulentMixing`, `noMixing`). Their source terms come from run-time-selectable sub-models under
  `populationBalanceSubModels/` (aggregation & coalescence kernels, breakup kernels, daughter
  distributions, growth, nucleation, diffusion, collision, environment mixing) and `mixingSubModels/`.
- **`realizableOdeSolver`** — adaptive ODE integration of the moment source terms that preserves
  moment realizability; configured via the `odeCoeffs` sub-dictionary.

`mappedList` / `mappedPtrList` (header-only, included via `-I../../mappedList`) are the indexing
backbone: they give named access to moments/nodes by moment order (`(0)`, `(1 0)`, `(2 1 0)`, ...)
with O(1) access by index and O(log n) by key. Moment orders and node indexes are declared as
`labelListList`s in `quadratureProperties`, and multi-dimensional orders are packed into a single
label key — bugs in this packing/unpacking are a recurring source of subtle breakage.

Nearly every model layer uses OpenFOAM's `declareRunTimeSelectionTable`/`addToRunTimeSelectionTable`
pattern. Adding a model means: new class in the right sub-directory, `addToRunTimeSelectionTable`,
add the `.C` to the library's `Make/files`, and it becomes selectable by keyword in the case
dictionary.

## Case setup conventions

- `constant/quadratureProperties.<quadratureName>` (or plain `quadratureProperties` when unnamed)
  selects `fieldMomentInversion`, the `basicMomentInversion`/`extendedMomentInversion` sub-dictionary,
  a `momentAdvection` sub-dictionary (which carries its *own* `basicMomentInversion` used for face
  reconstruction), and the `moments (...)` / `nodes (...)` lists of orders/indexes.
- `constant/populationBalanceProperties` selects `populationBalanceModel` plus its `<model>Coeffs`
  dictionary with the on/off switches and kernel sub-dictionaries. Similarly
  `constant/mixingProperties` for mixing.
- Moment fields are named `moment.<order>.<quadratureName>`, e.g. `moment.3.populationBalance`,
  and live in `0.orig/`; `Allrun` scripts copy `0.orig` → `0` before running.
- Multi-stage tutorials swap `system/controlDict.flow` / `controlDict.pbe` into `system/controlDict`
  between stages (solve the flow field first, then freeze it and solve the population balance).

## Applications

`applications/solvers` is grouped by physics: `populationBalance` (pbeFoam, pbeTransportFoam,
buoyantPbePimpleFoam), `mixing` (mixingFoam, mixingTransportFoam), `compressible`
(explicitRhoFoam + compressiblePbeTransportFoam, sharing a `ButcherTable` library),
`velocityDistributionTransport` (vdfTransportFoam, diluteVdfTransportFoam,
oneWayCoupledVdfTransportFoam + their own phaseModel/interfacial/turbulence libraries), and
`multiphase` (polydisperseBubbleFoam, denseAGFoam + twoPhaseSystem). Several of these groups
carry *private* copies of interfacial-model and phase-turbulence libraries — a change to one
does not propagate to the other; check whether the sibling copy needs the same fix.

Solvers are thin: they build a `populationBalanceModel`/`mixingModel` in `createFields.H` and
call `->solve()` inside the PIMPLE loop. Utilities in `applications/utilities` handle case
preparation and post-processing (`generateMoments` builds moment fields from a distribution,
`computeMoments`, `reconstructPointDistribution`, `errorEstimator`, `phaseMeanVelocityForce`).

## Conventions

- OpenFOAM coding style: 4-space indent, 80-column limit, `Foam::` namespace, the standard banner
  comment, `//- ` doc comments, and the `* * *` separator comments. Match the surrounding file.
- Every file keeps the OpenQBMM/OpenFOAM GPL header with its existing copyright lines — preserve
  and extend them rather than replacing.
- Commit messages use OpenFOAM prefixes: `BUG:`, `ENH:`, `DEV:`, `COMP:`, `STY:`, `PRJ:`.
- Version-dependent code is guarded with `#if (OPENFOAM >= NNNN)` (e.g. `2206` for `-lthermoTools`,
  `2512` for API changes) in both `Make/options` and sources — keep older versions compiling.
