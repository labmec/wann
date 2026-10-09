# Two-well example mesh

From the `wann` project root:

```sh
python3 geoMW/run_all.py
# Explicit input (paths are relative to geoMW, as in geoNew):
python3 geoMW/run_all.py ../input/twoWell.json
# Original finer outer-reservoir sizing:
python3 geoMW/run_all.py --reservoir-size-scale 1
# Output override (relative to the current working directory):
python3 geoMW/run_all.py ../input/twoWell.json --output /tmp/twoWell.msh
```

Requires Python 3 and the `gmsh` executable (`--gmsh` selects another executable).
Uses `geoNew/run_all.py` for JSON parsing, parameter defaults and material IDs.
The default input is `input/twoWell.json` (singular in this checkout).
Without `--output`, writes `MeshData.file` relative to the JSON directory.
The supplied input now selects `input/twoWell.msh`.

This example supports exactly two straight wells in a `pill` or `box` reservoir.
Their axis heels are `(0, -W/4, 0)` and `(0, W/4, 0)`; each toe is at
`(well.length, heel.y, 0)`. Thus their separation is `W/2`. The old `height`
and `eccentricity` inputs are ignored for this temporary placement convention.
The 1D well elements and endpoint physical points remain on the cylinder rim,
following `geoNew`'s coupling convention; the **axes** have `z = 0`.

The pill follows geoNew's stadium shape: its straight section spans the longest
well, with end-cap radius `W/2`. As in geoNew, `ReservoirData.length` controls
mesh sizing, rather than the pill's total geometric length. The box uses the
specified reservoir length, centered on the longest well's midpoint.

Each well uses a translated copy of geoNew's structured near-well topology,
with its own six physical IDs. A reservoir surface with two rectangular holes
is recombined into quadrilaterals and extruded into three layers of
hexahedra. Delaunay meshing, Blossom recombination, topology optimization
and smoothing produce predominantly hexahedral cells. A positive recombination
quality threshold leaves a few prisms where forcing quadrilaterals would
create degenerate corners. The near-well regions contain
hexahedra and match those hole interfaces. All reservoir cells share
`ReservoirData.matid`; cap-rock and farfield IDs also come from the JSON.
Well IDs must be distinct across both wells and from reservoir IDs.

The outer reservoir defaults to a mesh size multiplier of 2 for initial tests.
Use `--reservoir-size-scale` to adjust it: larger values give coarser elements;
`1` restores the original sizing. This multiplies the minimum and maximum
background-field sizes, while preserving the near-well discretization and
interface node counts.

Generation uses temporary intermediate files. Intermediate ASCII Gmsh 2.2 meshes are combined into a mesh that
preserves physical tags, offsets intermediate geometric entity tags, and merges
coincident nodes at `1e-8` coordinate precision. The combined mesh is exported through Gmsh to **ASCII MSH 4.1**, as required
by NeoPZ's `TPZGmshReader`. This conversion preserves nodes, connectivity and
physical IDs without remeshing. The driver rejects Gmsh error
logs even when the executable returns zero. Before export, the driver checks all four vertical sides of both well boxes.
Hexahedron Jacobians must also be positive at all eight corners; this catches
quads with collinear corners that integration-point checks can miss.
Every interface face must share node IDs between exactly two volume cells,
one inside the well box and one outside; the face count must match the
structured discretization. The supplied example has 228 shared faces per well.
Existing output is replaced only
after successful generation and merging.

Validation of the supplied pill example: 8,840 nodes, 10,173 elements,
including 7,677 connected reservoir volume cells: 7,653 hexahedra and 24 prisms. The original sizing produced
21,632 total elements and 14,502 volume cells. All 15 physical groups are
present; every exterior volume face has its boundary tag; interface faces share
nodes; Gmsh Jacobians are positive. This validates mesh construction, without
running the reservoir simulator.
