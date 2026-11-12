# Paper targets and implementation limits

See the [repository README](../README.md) for the three experiment commands.

## Adaptivity

`TPZWannAdaptivityTools` builds with `WANN_lib`. The adaptive target
prescribes a uniform flux on each well wall. It divides the net endpoint rate
by the cylindrical area `2 * pi * radius * length` from the JSON. With the bundled
input, it writes the refinement log to the working directory.
Refinement propagates around each well segment using the well's own axis
and includes its wall faces, coupling elements, and adjacent reservoir cells.

- `PragerSynge()` compares the mixed and H1 reservoir fluxes and adds a
  balance term for zero source. It stores indicators by geometric index and
  returns their Euclidean norm. Passing zero threads runs one worker.
- `AdaptivityProcess()` marks positive indicators at or above
  `theta * maximum`, enforces well compatibility, smooths the marks,
  rechecks compatibility, and records and divides the marked elements. The
  paper target uses `theta = 0.5` and compares the total scaled error against
  `1e-4`. Physical error values are obtained by dividing by
  `sqrt(pressureScale)`.

## Assumptions to keep visible

The geometry verifier and reservoir VTK exporter assume a 3D reservoir. The
coupling groups straight wells by their own heel-to-toe axes. The bundled
Ozkan generator creates a single x-aligned well. The mesh builders use the first
configured fluid, `fluid[0]`, for the single-phase solves.

The detailed material/output, copy, residual, and line-search limitations are
listed in [materials](materials.md), [nonlinear wells](nonlinear-wells.md), and
[post-processing](postprocessing.md).
