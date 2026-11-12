# Paper targets and implementation limits

See the [repository README](../README.md) for the three experiment commands.

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
