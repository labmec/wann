# Legacy code and current limits

The repository contains active library code alongside research experiments and
older drivers. The following observations come from current source and CMake
files; they are not a claim that every executable was built or numerically
validated as part of documenting the project.

## Drivers

| CMake target | Source | Current role/status |
|---|---|---|
| `twoWells` | `targets/twoWells.cpp` | Uses current multi-well data and `CreateMixedCompMesh()`; reads JSON, transforms well BC signs, solves mixed model, restores signs, and exports VTK |
| `wann3dNL` | `targets/wann3dNL.cpp` | Historical mixed/H1 comparison; still uses single-object well fields, `CreateMultiphysicsCompMesh()`, old diagnostic signatures, and an adaptivity include |
| `radial` | `targets/radialFlow.cpp` | Radial exact-solution/refinement experiment with local mesh factories; still uses old single-well configuration access |
| `visualize` | `targets/visualizeMesh.cpp` | Near-well/outer-mesh inspection experiment; still references removed material-ID fields |

`targets/needRefactor/` contains older adaptive and goal-oriented drivers not listed
in the current target CMake file. Treat these files as historical examples requiring
migration before reuse. Listing a target in CMake does not establish compatibility
with the present `ProblemData` definition.

## Standalone checks

`tests/CMakeLists.txt` defines manually run executables, not registered CTest tests:

| Target | Purpose indicated by the source |
|---|---|
| `hanging-nodes-test` | Mixed/H1 comparison and refinement/hanging-node exploration |
| `test-divide` | Geometry division and refinement experiments |
| `test-sideOrient-1d` | One-dimensional side/flux orientation inspection |
| `test-nonlinear-3d` | Linear versus Newton Darcy assembly experiments |
| `testDarcyNL` | Darcy Newton comparison/check driver |
| `test-nonlinear-1d` | One-dimensional mixed/H1 nonlinear well experiments |

They should be assessed individually before interpreting them as automated
regression coverage. Some use historical APIs and contain old solver variants.

## Historical adaptivity class

`src/needRefactor/TPZWannAdaptivityTools` is excluded by the nonrecursive library
source glob. Its implementation still refers to removed fields such as `EDomain`,
`ECurveWell`, and `ESurfWellCyl`.

Its intended method responsibilities are documented in its header:

- `PragerSynge()` compares mixed and H1 reservoir fields and adds a balance term;
  local indicators are stored by geometric index and the returned aggregate is
  their Euclidean norm. The source is currently zero and permeability is diagonal.
- `GoalOriented()` accumulates signed primal/dual goal contributions while storing
  absolute local indicators. The well contribution is disabled in its dispatcher.
- `ErrorEstimation()` is an older scalar-permeability mixed/H1 comparison with a
  well special case and diagnostic log output.
- `SimpleMarking()` compares each indicator to `theta * maximum`.
  `BulkMarking()` sorts indicators and marks until their sum reaches
  `theta * total`; it uses values directly, not squared values.
- `AdaptivityProcess()` marks, applies well compatibility, smooths, rechecks
  compatibility, records geometric indices, and divides marked elements.
- `ForcingFunctionWellbore()` is a zero-valued placeholder.
  `GoalContributionInterface()` is declared without an implementation.

Estimator routines increment the requested thread count before launching workers,
so an argument of zero gives one worker. Mesh/solution references and material IDs
must correspond across the compared meshes. These historical routines do not
supply a supported current adaptivity workflow until their API is migrated.

## Assumptions to keep visible

Tensor support for 1D/2D input does not automatically generalize every simulation
path. The geometry verifier and reservoir VTK exporter assume three dimensions;
H1 coupling currently selects the first x-aligned well; mixed-well assembly uses
only the x component of flux. Current mesh builders are single-phase and use
`fluid[0]`, regardless of the number of fluids parsed.

The detailed material/output, copy, residual, and line-search limitations are
listed in [materials](materials.md), [nonlinear wells](nonlinear-wells.md), and
[post-processing](postprocessing.md). These documents intentionally describe the
current implementation without modifying those paths.
