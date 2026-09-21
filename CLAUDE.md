# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

OpenQBMM is an OpenFOAM **module** (not a standalone project) implementing Quadrature-Based
Moment Methods (QBMM) for polydisperse multiphase flows: population balance, turbulent mixing,
and velocity-distribution (kinetic) transport. It targets the **openfoam.com (ESI) line**;
`master` and `development-openfoam.com` are the live branches. CI (`.github/workflows/main.yml`)
builds the libraries, the applications and the unit tests against OpenFOAM 2412 and runs the unit
tests.

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
`test/Allwmake` deliberately skips `conditionalMomentInversion` (commented out, needs updating for
the current structure).

## Cases, validation and regression

`tutorials/` and `validation/` are per-solver OpenFOAM cases, each with `Allrun`/`Allclean`.
Allrun scripts source `AllrunFunctions` from the top of the sources and use `runStep` /
`runStepParallel` (which stop the case on a failed or skipped step) and `plotCase`; parent
scripts run their cases with `runCase` and end with `reportCases`.

- **Regression.** A case with `regression/baseline` is checked by `bin/qbmmRegression check <case>`
  against a summary of the fields of its last time (min/max/mean/rms per field, from
  `bin/qbmmFieldSummary`), relative to the span, at tolerance 1e-8 unless
  `regression/tolerance` says otherwise. `qbmmRegression record <case>` re-records; that is a
  decision to take case by case, never done silently. 60 cases carry baselines.
- **Reference data.** A validation case declares its figures and its comparisons with reference
  data in one OpenFOAM-style dictionary, `validation/validationDict` (format in
  `bin/python/qbmm/validation.py`): each curve once, in its figure, and a computed curve names in
  `compare` the reference curves it is held against with a tolerance on the rms difference (as a
  fraction of the reference's scale) or `report` (measured, not judged).
  `bin/qbmmCompareReference <case>` measures them (pure Python, no dependencies); `bin/qbmmPlot
  <case>` draws the figures with matplotlib into `validation/figures/` (SVG + PDF), and is skipped
  with a message where matplotlib is missing (`QBMM_PYTHON` names another interpreter). Anything
  computed from logs or probes before plotting is done in the Allrun or a case script, not in the
  dictionary. OpenFOAM `raw` sets write their fields sorted by name, not in the order of
  `sampleDict`, so column numbers follow the name of the written file.
- **Tiers.** `bin/qbmmValidate -quick|-baselines|-complete` runs a tier of cases plus the unit
  tests, in copies under `$FOAM_RUN/qbmmValidate/`, as many at once as there are cores (a parallel
  case counts its subdomains), and writes a report. Quick is ~2 min, baselines and complete
  ~2 h on 22 cores. `validation/validationTiers` lists the `long` cases (skip quick) and the
  `manual` ones (never run: polydisperseBubbleFoam, the Taylor-Couette cases); a new short case
  joins quick by itself. `-list` prints a tier.
- **Report.** `bin/qbmmReport -site <GH Pages export> [-tests <test/Allrun log>]` collects the
  state the cases were left in (nothing is run) into `<site>/validation/`: the dashboard page
  (header and footer taken from the site's `index.html`, stylesheet `../industry.css`),
  `results.json`, a per-commit `history.json` and SVG figures; it requires matplotlib.
  `-output <dir> -style <css>` writes a standalone copy, `-refresh` re-renders the page from the
  results it holds. Figures use the palette of the OpenQBMM site (industry design system).

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
  - `multivariate`: CHyQMOM, CHyQMOM+, sizeCHyQMOM, TensorProduct, conditionalMeans (also
    selected as `monoKinetic`: a primary coordinate with the means of velocities or scalars
    conditioned on it, supports taken from the moment set, scalars clipped to theirs; as many
    nodes, and first-order moments of each conditioned coordinate, as the univariate rule of
    its `basicQuadrature` builds, GQMOM included, with a warning beyond 2N - 1 nodes),
    conditional.
  These operate on plain `momentSet`s — no mesh, no fields. This is what the unit tests exercise.
- **`fieldMomentInversion`** — applies a `momentInversion` cell-by-cell over an `fvMesh`
  (`basicFieldMomentInversion`, `basicMultivariateFieldMomentInversion`, formerly
  `basicVelocityFieldMomentInversion`, still read with a warning). This is the bridge from
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
  boundary handling, plus `firstOrderFluidVelocity`, first-order upwind for a quadrature whose nodes all
  move with the fluid, its flux named by `phi` in the advection dictionary).
- **`PDFTransportModels`** — abstract base tying advection + source terms together
  (`univariatePDFTransportModel`, `velocityPDFTransportModel`).
- **`populationBalanceModels`** / **`mixingModels`** — the concrete physics on top of PDF transport,
  selected at run time (`univariate`, `velocity`, `sizeVelocity`, `mixingPopulationBalance`, `none`;
  and `turbulentMixing`, `conditionalMixing` — a mixture fraction with the compositions
  conditioned on it, both on [0,1], carried by the fluid and mixed by the same kernels —,
  `noMixing`). Their source terms come from run-time-selectable sub-models under
  `populationBalanceSubModels/` (aggregation & coalescence kernels, breakup kernels, daughter
  distributions, growth, nucleation, diffusion, collision, environment mixing) and `mixingSubModels/`.
- **`realizableOdeSolver`** — adaptive ODE integration of the moment source terms that preserves
  moment realizability; configured via the `odeCoeffs` sub-dictionary.

`mappedList` / `mappedPtrList` (header-only, included via `-I../../mappedList`) are the indexing
backbone: they give named access to moments/nodes by moment order (`(0)`, `(1 0)`, `(2 1 0)`, ...)
with O(1) access by index and by order. The index is `mappedListKey`: the orders of every element
in a flat list and an open-addressed table on an exact key (the orders packed into 64 bits, or a
hash checked against the stored orders when they do not fit), so two different orders are never
confused. Orders followed by zeros are the same element (`(1)` is `(1 0 0)`); `found()` is strict
on the number of dimensions. `mappedListKey::ordersToWord` gives the word used in field names.

The moment set of a quadrature is built by its inversion, not listed in the case:
`fieldMomentInversion::momentSet` (and `momentSupports`) give the orders, node indexes and
supports from the few counts the inversion dictionary carries. Conditioned inversions
(sizeCHyQMOM, conditionalMeans, conditional) share `Vandermonde::conditionalInverse`.

Nearly every model layer uses OpenFOAM's `declareRunTimeSelectionTable`/`addToRunTimeSelectionTable`
pattern. Adding a model means: new class in the right sub-directory, `addToRunTimeSelectionTable`,
add the `.C` to the library's `Make/files`, and it becomes selectable by keyword in the case
dictionary.

## Case setup conventions

- `constant/quadratureProperties.<quadratureName>` (or plain `quadratureProperties` when unnamed)
  selects `fieldMomentInversion` and its inversion sub-dictionary, which carries the counts the
  moment set is built from: `basicMomentInversion { univariateMomentInversion ...; nMoments N; }`
  (plus `nAdditionalQuadraturePoints` for GQMOM; Radau and Lobatto add one and two nodes), or
  `basicMultivariateMomentInversion { type CHyQMOM|CHyQMOMPlus|sizeCHyQMOM|monoKinetic;
  nVelocityDimensions D; }` (plus `nSizeMoments` for sizeCHyQMOM and monoKinetic, and
  `transportedMoments` when the case transports another closure's moments). The old names
  `basicVelocityFieldMomentInversion` / `basicVelocityMomentInversion` are read with a warning;
  `bin/qbmmConvertQuadratureProperties` renames them. The supports of a `conditionalMeans`
  quadrature are the model's, not the dictionary's. A `momentAdvection` sub-dictionary carries
  its *own* inversion used for face reconstruction, and `timeIntegration Euler|SSPRK2` (default
  Euler; SSPRK2 needs `ddtSchemes Euler` and no under-relaxation of the moment equations, which
  the transport models check).
  Cases that still list `moments (...)` / `nodes (...)` are read with a warning;
  `bin/qbmmConvertQuadratureProperties <case>` converts them.
- `constant/populationBalanceProperties` selects `populationBalanceModel` plus its `<model>Coeffs`
  dictionary with the on/off switches and kernel sub-dictionaries. Similarly
  `constant/mixingProperties` for mixing.
- Moment fields are named `moment.<order>.<quadratureName>`, e.g. `moment.3.populationBalance`;
  a moment of several dimensions joins its orders with underscores, `moment.2_1_0.particles`, and
  so do node names (`weight.node0_1_2.particles`). Cases written with the old concatenated names
  (`moment.210.particles`) are still read with a warning; `bin/qbmmConvertMomentNames <case>`
  converts them. Both converters are temporary. Moment fields live in `0.orig/`; `Allrun`
  scripts copy `0.orig` → `0` before running.
- Multi-stage tutorials swap `system/controlDict.flow` / `controlDict.pbe` into `system/controlDict`
  between stages (solve the flow field first, then freeze it and solve the population balance).

## Applications

`applications/solvers` is grouped by physics: `populationBalance` (pbeFoam, pbeTransportFoam,
buoyantPbePimpleFoam), `mixing` (mixingFoam, mixingReactingFoam), `compressible`
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
