# Export and well diagnostics

`TPZWannPostProcTools` operates on solved meshes. Geometric element references can
point to different computational meshes when mixed and H1 meshes share geometry;
methods that evaluate through these references reload them before use where
implemented.

## VTK output

`WriteWellboreVTK()` exports `Pressure` and `Flux` for each well material. Mixed
well output also includes `Divergence`. Filenames include the configured prefix, `Hdiv` or `H1`, and
the zero-based well index. The routine identifies a mixed solve by dynamically
casting the mesh to `TPZMultiphysicsCompMesh`.

`WriteReservoirVTK()` exports `Pressure` and `Flux` using the configured reservoir
prefix and formulation suffix. It selects dimension three and passes the
`ToCylindrical` option as the VTK safety-check flag. It is a coupled 3D reservoir
exporter, rather than a generic dimension-independent material-output function.
Both exporters use the configured VTK resolution and post-processing thread count.

`WriteVTKs()` calls both exporters.

## Well profiles for the Ozkan comparison

`WriteWellProfilesCSV()` samples the mixed well solution along the normalized
heel-to-toe coordinate. For the bundled Ozkan input it writes three files with the
columns `x,Curve1`, where `x = s/L`.

| File | Value in `Curve1` |
|---|---|
| `fluxOzkanWell.csv` | Flow magnitude in STB/day |
| `divergenceOzkanWell.csv` | Dimensionless inflow magnitude, `L * abs(dQ/ds) / abs(Q_heel)` |
| `pressureOzkanWell.csv` | Far-field minus well pressure in psi |

`wann3dNL` uses four subdivisions per axial element. At interior interfaces it
samples the toe-side element. The inflow column is empty when the heel rate is zero.

## Segmented reservoir inflow

`ComputeWellFluxes(mesh, data, wellIndex, nSegments)` returns one integrated signed
normal flux for each equal-length axial bin of the selected well. It validates
the well index, positive segment count, distinct heel/toe BC IDs, and finite
nonzero endpoint distance.

It locates reservoir faces neighboring the selected well surface and projects each
face centroid onto the heel-to-toe axis. The entire face contribution is assigned
to the centroid's bin; it is not subdivided where a bin boundary crosses the face.
It evaluates the reservoir material’s `Flux` at the corresponding point in the
neighboring volume element. It evaluates the wall Jacobian and the outward
reservoir normal at each quadrature point, so positive flux is into the well.

## Productivity index

`ProductivityIndex()` looks for endpoint-adjacent well elements, reads their flow
and pressure, and compares pressure to reservoir BC `surface_farfield`. It unscales
the pressure difference according to

```text
deltaP_Pa = (endpoint pressure - farfield pressure) / pressureScale
PI        = endpoint flow / deltaP_Pa
```

It converts the result from m³/(s Pa) to STB/(day psi) using the existing factors
6.28981 barrels/m³, 86400 seconds/day, and 6894.76 Pa/psi. The pressure difference
and flow remain signed; the reference-pressure offset cancels.

The routine evaluates the endpoint adjacent to the heel on either local side of
the well element. It assumes a nonzero pressure difference.
