# Export and well diagnostics

`TPZWannPostProcTools` operates on solved meshes. Geometric element references can
point to different computational meshes when mixed and H1 meshes share geometry;
methods that evaluate through these references reload them before use where
implemented.

## VTK output

`WriteWellboreVTK()` exports `Pressure`, `Flux`, `Divergence`, and `GradPressure` for
each well material. Filenames include the configured prefix, `Hdiv` or `H1`, and
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
heel-to-toe coordinate. It exports flow magnitude in STB/day, dimensionless inflow
magnitude, and far-field pressure drop in psi.

## Segmented reservoir inflow

`ComputeWellFluxes(mesh, data, wellIndex, nSegments)` returns one integrated signed
normal flux for each equal-length axial bin of the selected well. It validates
the well index, positive segment count, distinct heel/toe BC IDs, and finite
nonzero endpoint distance.

It locates reservoir faces neighboring the selected well surface and projects each
face centroid onto the heel-to-toe axis. The entire face contribution is assigned
to the centroid's bin; it is not subdivided where a bin boundary crosses the face.
It obtains reservoir flux from atomic mixed mesh zero or the H1 material's flux
output, and uses the geometric face normal for the sign. Opposite normal
orientation reverses the signed result.

The implementation constructs its face Jacobian/normal before calling
`intrule->Point()` to load that iteration's quadrature coordinates. It should not
be assumed to evaluate curved-face geometry correctly at every quadrature point.
This existing ordering is left unchanged.

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

The current loop accepts a heel neighbor at local side zero or a toe neighbor at
local side one, then overwrites the stored flow and pressure for every matching
element. The returned value is therefore not guaranteed to be a heel-only value
and can depend on element order. There is no explicit zero-pressure-difference
check. These limitations should be considered when interpreting the diagnostic.
