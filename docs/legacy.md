# Legacy code and current limits

The repository contains active library code alongside research experiments and
older drivers. The following observations come from current source and CMake
files; they are not a claim that every executable was built or numerically
validated as part of documenting the project.

## Drivers

| CMake target | Source | Current role/status |
|---|---|---|
| `twoWells` | `targets/twoWells.cpp` | Uses current multi-well data and `CreateMixedCompMesh()`; reads JSON, transforms well BC signs, solves mixed model, restores signs, and exports VTK |
| `wann3dNL` | `targets/wann3dNL.cpp` | Nonlinear mixed and H1 coupled solves with per-well flux and productivity. Default input is `ozkan1999.json`. |
| `radialFlowTest` or `radial` | `targets/radialFlow.cpp` | Radial convergence for mixed degree 1 and H1 degrees 1 and 3 with the default input `radialTest.json`. Runs five mesh levels by default. Pass `radialTest 1` for a two-level check. |
| `adaptiveRefinement` | `targets/adaptiveRefinement.cpp` | Mixed/H1 reservoir error estimation with well-compatible refinement. Default input is `ozkan1999adaptive.json`. |
| `visualize` | `targets/visualizeMesh.cpp` | Exports the bundled near-well and outer mesh parts. Pass a JSON input to export its configured mesh. |

`targets/needRefactor/` contains older research drivers not listed
in the current target CMake file. Treat these files as historical examples requiring
migration before reuse. Listing a target in CMake does not establish compatibility
with the present `ProblemData` definition.

## Standalone checks

`tests/CMakeLists.txt` defines manually run executables, not registered CTest tests:

| Target | Purpose indicated by the source |
|---|---|
| `hanging-nodes-test` | Mixed/H1 comparison and refinement/hanging-node exploration |
| `test-divide` | Geometry division and refinement experiments |
| `test-sideOrient-1d` | Verifies the signed heel-to-toe flux in a one-dimensional mixed Darcy solve |
| `test-nonlinear-3d` | Mixed and H1 Newton solves for a unit-cube Darcy problem |
| `testDarcyNL` | Mixed and H1 Newton solves on the historical reservoir mesh |
| `test-nonlinear-1d` | One-dimensional mixed/H1 nonlinear well experiments |

They should be assessed individually before interpreting them as automated
regression coverage. Some use historical APIs and contain old solver variants.

## Adaptivity

`TPZWannAdaptivityTools` builds with `WANN_lib`. The adaptive target
prescribes a uniform flux on each well wall from its net endpoint rate and
wall area. Its refinement log is written to the working directory.
Refinement propagates around each well segment using the well's own axis
and includes its wall faces, coupling elements, and adjacent reservoir cells.

- `PragerSynge()` compares the mixed and H1 reservoir fluxes and adds a
  balance term for zero source. It stores indicators by geometric index and
  returns their Euclidean norm. Passing zero threads runs one worker.
- `AdaptivityProcess()` marks positive indicators at or above
  `theta * maximum`, enforces well compatibility, smooths the marks,
  rechecks compatibility, and records and divides the marked elements.

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
