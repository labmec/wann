# Mesh generation

Run the driver with the simulation JSON:

```sh
python3 geoTest/run_all.py ../input/ozkan1999.json
```

Relative JSON paths are resolved against `geoTest`, regardless of the current
working directory. With no argument, the driver uses `input/ozkan1999.json`.
The final mesh is written to `MeshData.file`, relative to the JSON directory.

The driver accepts `WellboreData` as either an object or a one-element array.
The current geometry templates construct one straight well at the requested rotation angle;
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
All four formats (`box`, `pill`, `ball`, `nearWellbore`) use these parameters.

Older JSON files without the required material IDs must be updated before use.
Direct execution of the parameterized `.geo` files requires supplying the tag
parameters; missing tags produce an error. `reservoirBoxStruct.geo` is a standalone
mesh experiment and is not called by the driver.

Control the rotation with `--angle`, in degrees from +x toward +y:

```sh
python3 geoTest/run_all.py ../input/ozkan1999.json --angle 0
python3 geoTest/run_all.py ../input/ozkan1999.json --angle 45
python3 geoTest/run_all.py ../input/ozkan1999.json --angle 90
```

The default remains 90 degrees. The well and surrounding reservoir rotate
about the z-axis through the origin; height and dimensions are preserved.
All templates use `rotation.geo`, with the mapping
`(x, y, z) -> (x*cos(a) - y*sin(a), x*sin(a) + y*cos(a), z)`.
Side comments describe the local frame before rotation.
For direct Gmsh execution, supply `-setnumber well_angle 45` along with
physical tags. The standalone structured experiment accepts this parameter too.

As in `geoNew`, the output path comes from `MeshData.file`. Use a separate
JSON output path if you want to keep both x- and y-aligned meshes.
