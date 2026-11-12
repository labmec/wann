# Mesh generation

## Ozkan mesh

Run the driver with the simulation JSON:

```sh
python3 meshing/run_all.py ../input/ozkan1999.json
```

Relative JSON paths are resolved against `meshing`, regardless of the current
working directory. With no argument, the driver uses `input/ozkan1999.json`.
The final mesh is written to `MeshData.file`, relative to the JSON directory.

The driver accepts `WellboreData` as either an object or a one-element array.
The current geometry templates construct one straight well along the x-axis;
multiple-well arrays are rejected rather than silently generating only one well.

Physical tags are supplied through `-setnumber` to the shared `params.geo`:

| Physical group | JSON source |
| --- | --- |
| `curve_wellbore` | `WellboreData.matid` |
| `surface_wellbore_cylinder` | `WellboreData.matidSurf` |
| `surface_wellbore_toe` | `WellboreData.matidToeSurf` |
| `surface_wellbore_heel` | `WellboreData.matidHeelSurf` |
| `point_heel`, `point_toe` | Matching entries in `WellboreData.BCs`, field `matid` |
| `volume_reservoir` | `ReservoirData.matid` |
| `surface_farfield`, `surface_cap_rock` | Matching entries in `ReservoirData.BCs`, field `matid` |

The required IDs must be distinct positive integers. Both intermediate meshes
receive the same reservoir and cap-rock IDs, so merging preserves these groups.
The driver generates the pill-shaped reservoir used by the Ozkan experiments.

Regenerate the adaptive refinement history whenever you generate a new Ozkan
mesh. The history refers to geometric element indices in the original mesh.

## Radial mesh

From the repository root, run the radial Gmsh script from the `meshing` directory.

```sh
(cd meshing && gmsh radialTest.geo -)
```

The script writes `input/radialTest.msh`.
