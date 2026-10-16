# Configuration, mesh construction, and solving

## Configuration

Construct `ProblemData`, then call `ReadJson(filename)` before using its fields.
The filename is resolved relative to the CMake-configured `INPUTDIR`. The parser
accepts comments in JSON. It checks required fields with `DebugStop()`; malformed
JSON or conversions can also throw through the JSON library.

| Stored structure | Content and parsing behavior |
|---|---|
| `BoundaryData` | Material ID, integer type, and value; name is the enclosing map key |
| `DomainData` | Name, material ID, permeability vector, and named boundaries |
| `WellboreData` | Domain metadata, cylinder/toe/heel surface IDs, radius, length, three-entry eccentricity, and BCs; accepts one object or an array |
| `ReservoirData` | Domain metadata; scalar permeability replicated into three diagonal entries or an array of exactly three; geometry dimensions and porosity |
| `FluidData` | Name, viscosity, density; accepts one object or an array |
| `MeshData` | Mesh filename, custom/uniform/directional refinement requests, and cylindrical-map flag |
| `PostProcData` | VTK/training names and resolutions; optional thread count/verbosity default to zero |
| `NumericsData` | Threads, iteration/tolerance settings, orders, and pressure scaling; the object itself is required |

The constructor only sizes the reservoir permeability vector and reserves boundary
map capacity. It does not populate every scalar field. The inherited well
permeability vector is not filled by the current JSON reader. The current
approximation builders use `m_Fluid[0]`; reading an array does not implement a
multiphase constitutive model.

Default numerical settings used by `ReadJson()` are:

| Field | Default |
|---|---:|
| `nthreads` | 0 |
| `maxIterations` | 10 |
| `res_tol`, `corr_tol` | 1e-6 |
| `reservoirPorder` | 1 |
| `wellPorder` | 2 |
| `refPressure` | 1e5 |
| `pressureScale` | 1e5 |

After parsing, the reader calls three mutating methods:

1. `UpdateAuxiliaryMaterialIds()` finds the largest configured domain/surface/BC
   ID and assigns the next three IDs to the pressure skin, pressure interface,
   and H(div) boundary interface.
2. `UpdateReservoirBCsForWellbores()` inserts zero-Neumann reservoir boundaries
   on every well toe and heel surface. Their names are `no_flux_toe_<well name>`
   and `no_flux_heel_<well name>`.
3. `ApplyPressureScaling()` multiplies fluid viscosity and density by
   `pressureScale`; type-0 pressure BCs become
   `(physical pressure - refPressure) * pressureScale`. Other BC values are
   unchanged. This method is not idempotent.

The stored values after `ReadJson()` are therefore not all the original physical
input values. `Print()` reports the stored configuration. Do not apply pressure
scaling again to data already read/scaled.

## Geometry

`CreateGeoMesh()` performs the following operations:

1. `ReadMeshFromGmsh()` imports `m_Mesh.file` relative to `INPUTDIR`, then replaces
   codimension-one Gmsh boundary elements with `TPZGeoElBC` elements to establish
   consistent orientation.
2. `VerifyMesh()` checks 3D reservoir IDs and face neighbor counts, and 2D boundary
   neighbors. It returns **true for errors**, assigns diagnostic IDs 1001–1005,
   and causes the factory to write `gmeshWithErrors.vtk` and stop if errors exist.
3. `CreateCouplingEls()` places three auxiliary geometric elements on each well
   surface: pressure skin, pressure interface, and H(div) boundary interface.
4. `OrderWellIds()` finds the named `point_heel`/`point_toe` boundaries of each well
   and calls `OrderIdsSingleWell()` to assign increasing node IDs by axial group.

The factory currently does not execute cylindrical conversion or refinement from
`MeshData` flags. `ModifyGeometricMeshToCylWell()` is a separate x-aligned cylinder
helper; it projects nodes and converts selected maps and adjoining GeoBlend
geometry. `hRefinement()` divides a list of leaf geometric indices.
`RefineFromFile()` replays lists, one level per line, in the form
`count index0 index1 ...`.

`ComputeAxialCoordinate(point, origin, axis)` projects onto the normalized axis.
It requires finite nonzero axis length and three-component vectors. The coordinate
set helpers group nearly equal scalars with a caller-supplied tolerance; their
historical `X` names also apply to projected axial coordinates. `FindClosestX()`
is a tolerance lookup, not an unrestricted nearest-neighbor query.

## Mixed approximation

`CreateMixedCompMesh()` constructs the coupled approximation in this order:

1. Configure `TPZHDivApproxCreator` with the reservoir order and condensation off.
2. Insert the Newton mixed reservoir material and its BCs. Supply a diagonal
   tensor with entries `reservoir.perm[i] / fluid[0].viscosity`.
3. Insert null coupling materials and one physical mixed well material per well.
   The builders currently pass zero reference pressure and pseudo-resistivity;
   reservoir/well exchange is represented by the coupling spaces/interfaces.
4. Create atomic meshes. Index zero contains flux and index one pressure.
5. `AddPressureSkinElements()` creates disconnected quadrilateral skins, removes
   circumferential/face interior modes, and shares axial pressure connects around
   the circumference within each segment of each well. It uses the actual
   heel-to-toe projection, requires positive axial extent and `wellPorder >= 1`,
   and checks each skin has one well owner.
6. `AddWellboreElements()` adds disconnected line-pressure elements and H(div)
   line-flow/end-point elements for every well.
7. `EqualizeL2PressureConnects()` makes well pressure reuse the neighboring skin
   edge and nodal connects, respecting node-order reversal.
8. `AddHDivBoundInterfaceElements()` creates normal-flux trace elements, and the
   creator assembles the multiphysics mesh.
9. `AddInterfaceElements()` links skin-pressure and flux-trace elements with
   `TPZMultiphysicsInterfaceElement` and a Lagrange multiplier material.

`ComputeTotalMobility()` is a private declaration without an implementation or
call site. Mobility is currently supplied directly as `K/mu`, not through that
helper. `AddInterfaceElements()` does not use its `laglevel` argument.

## H1 approximation

`CreateH1CompMesh()` creates continuous pressure spaces with the Newton H1 reservoir
and a well material for each configured well. It optionally attaches a `TLaplaceExample1` exact solution
and forcing function to the well. It then updates line-interior order and calls
`EqualizeH1PressureConnects()`.

That helper identifies each skin's well through the neighboring cylinder material
ID and projects coordinates onto that well's straight heel-to-toe axis. It groups
vertex and axial-edge connects separately within each well, disables
circumferential/face modes, applies identity dependencies within coordinate
groups, and rebuilds solution blocks. Wells can have arbitrary straight axes;
sharing never crosses well ownership. Curved wells require a path coordinate
rather than this endpoint projection.

## Analysis lifecycle

For a compatible mixed mesh, the library sequence is:

```cpp
ProblemData data;
data.ReadJson("case.json");
TPZGeoMesh *geometry = TPZWannGeometryTools::CreateGeoMesh(&data);
TPZMultiphysicsCompMesh *mesh =
    TPZWannApproxTools::CreateMixedCompMesh(geometry, &data);

{
    TPZWannAnalysis analysis(mesh);
    analysis.SetProblemData(&data);
    analysis.Initialize();
    analysis.NewtonIteration();
    TPZWannPostProcTools::WriteVTKs(mesh, &data);
} // Analysis lifetime ends before its mesh is deleted.

delete mesh;
delete geometry;
```

This is a library lifecycle example, not a complete reproduction of every driver:
`twoWells.cpp` also applies and reverts its own well-BC sign transformation. Choose
BC signs for the material and geometric orientation being used.

`Initialize()` chooses threaded MKL sparse assembly when enabled, or skyline
assembly otherwise. It gathers type-1 BC IDs for mixed meshes or type-0 IDs for H1,
initializes their solution values, and excludes their connect equations from the
correction solve. H1 initialization assigns corner values; mixed initialization
includes a boundary-measure/ncorner factor. It selects an LDLt direct solver for
mixed systems and Cholesky for H1 systems.

`NewtonIteration()` is the entire iteration loop, despite its singular name. It
assembles `J` and `-R`, solves for a correction, applies the selected scale, and
loads the new solution into the computational mesh. Mixed updates also transfer
the multiphysics solution to atomic meshes. After the first step, convergence
requires both residual and scaled correction norms to satisfy their tolerances.
The function prints convergence/nonconvergence rather than returning a status.
`GetNumberOfIterations()` exposes the loop counter.

See [well and line-search details](nonlinear-wells.md) for diagnostic side effects
and regime differences.
