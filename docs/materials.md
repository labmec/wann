# Darcy materials and permeability

## Permeability provider

`TPZAnisotropicPermeability` has no space-dimension member. The Darcy material owns
its `fDim`; the permeability provider always evaluates tensors in three-dimensional
global-coordinate storage.

| Input | Returned tensor | Frame convention |
|---|---|---|
| Scalar `k` | `diag(k, k, k)` | Isotropic in global coordinates |
| 1 × 1 matrix `[k]` | `diag(k, k, k)` | Physical 1D input associated with global x |
| 2 × 2 matrix `K2` | `diag-block(K2, sqrt(det(K2)))` | Physical 2D input associated with global xy; zero z coupling |
| 3 × 3 matrix | The supplied tensor | Global coordinates |

The same normalization applies to constant tensors and tensors returned by an
anisotropic callback. `GetPermeability(coord, matrix)` resizes/replaces its output
with a 3 × 3 matrix. `GetInversePermeability()` returns the cached constant inverse
or evaluates and inverts a spatial tensor each time.

The scalar getter has a different contract: it returns a value only for
scalar-configured permeability. Matrix setters mark the field anisotropic, even
if their matrix happens to be isotropic. Querying such a field through the scalar
getter prints a diagnostic and returns zero; it does not return the determinant
average stored internally by the constant-tensor setter.

Call a setter before using the provider. Its constructor initializes matrix
members using their scalar constructor; it does not perform the explicit 3 × 3
allocation and inverse setup done by `SetConstantPermeability()`.

For a constant tensor, the class computes its determinant and caches its LU
inverse. Shape checks and selected positivity/determinant checks are not a full
symmetric-positive-definite validation. Spatial 3 × 3 callbacks also do not pass
through the constant setter's determinant check. Supply valid physical tensors
and positive scalar permeabilities. The `Read()` and `Write()` hooks currently
serialize no permeability state.

## Why the embedding preserves mixed multiplication

For a 2D xy problem the physical flux basis has columns `[vx, vy, 0]^T`. With

```text
K3 = [ K2   0 ]      inverse(K3) = [ inverse(K2)          0 ]
     [  0  kN ]                    [      0           1/kN ]
```

and positive `kN`, the out-of-plane block is invertible but contributes nothing
to `v_i^T inverse(K3) v_j`. The 1D case follows similarly for `[vx, 0, 0]^T`.
This is why the Darcy mixed material can retain its 3 × 3 multiplication.

For a tilted plane, express the tensor globally. If `T` contains two orthonormal
tangent columns and `n` is a unit normal, a suitable extension is
`K3 = T K2 T^T + kN n n^T`. A reduced tensor alone does not encode this orientation.
The 1D isotropic extension `k I3` acts as `k` along any tangent, but does not infer
or encode a separate directional anisotropy.

Example using the current API:

```cpp
TPZFMatrix<STATE> K2(2, 2, 0.);
K2(0, 0) = 4.;
K2(0, 1) = K2(1, 0) = 1.;
K2(1, 1) = 3.;
material.SetConstantPermeability(K2);

TPZFMatrix<STATE> inverse;
material.GetInversePermeability(x, inverse); // 3 x 3
```

The coupled mesh builders currently supply diagonal `reservoir.perm / viscosity`.
Thus their tensor acts as hydraulic mobility, even though the provider calls it
permeability. The provider itself applies no viscosity factor.

## Shape functions and local blocks

For `TPZMixedDarcyAnisotropic` and `TPZWannMixedDarcyNL`, the primary data are:

| Expression | Size at a volume integration point | Meaning |
|---|---|---|
| `datavec[0].fDeformedDirections` | 3 × nFlux | Complete physical flux basis vectors, one per column |
| `datavec[0].divphi` | nFlux × 1 | Scalar divergence of each flux basis vector |
| `datavec[1].phi` | nPressure × 1 | Scalar pressure shape values, one per row |
| `datavec[0].sol[0]` | Physical flux vector | Current flux at this point |
| `datavec[0].divsol[0][0]` | Scalar | Current flux divergence |
| `datavec[1].sol[0][0]` | Scalar | Current pressure |

Pressure shape values depend on the integration point, element topology, and
polynomial space. The dimension does not introduce additional pressure components.
Boundary assembly instead uses `datavec[0].phi` as the normal-flux trace basis;
it is not interchangeable with the volume vector basis.

With local unknowns ordered `[q coefficients, p coefficients]`, the mixed volume
matrix is assembled from `A_ij = v_i^T K^{-1} v_j` and the two negative divergence
coupling blocks. The pressure-pressure block is zero for these Darcy materials.
NeoPZ supplies the physical integration weight; the routines add contributions
to existing `ek`/`ef` storage rather than resetting it.

The H1 Darcy material uses scalar `data.phi` and multiplies its 3 × 3 tensor directly
by `data.dphix`. This requires a compatible derivative matrix and coordinate
frame. Unlike its flux post-processing, volume assembly does not explicitly call
`Axes2XYZ`. Tensor embedding alone does not establish that every lower-dimensional
H1 element or rotated frame is compatible with this direct multiplication.

## Linear and Newton assembly

| Material | Volume right-hand side | Boundary behavior |
|---|---|---|
| `TPZDarcyAnisotropic` | Source load | Type 0 pressure penalty; type 1 load; type 2 Robin |
| `TPZMixedDarcyAnisotropic` | Source in pressure rows | Type 0 pressure load; type 1 normal-flux penalty; types 2/4/5 special terms |
| `TPZWannDarcyNL` | Negative pressure residual plus source | Type 0 handled by analysis; type 1 negative load; other types rejected |
| `TPZWannMixedDarcyNL` | Negative flux/pressure residual; no current source term | Type 0 pressure load; type 1 handled by analysis; no assembled Robin branch |

The `NL` names identify Newton assembly, not a pressure-dependent permeability
law. Their explicit permeability methods forward to the provider. Their stream
methods delegate to their base; they do not add serialization of the provider's
callbacks or cached tensor state.

## Post-processing and errors

Use each material's own `VariableIndex()` and `NSolutionVariables()`. Numeric IDs
are not interchangeable between H1, mixed Darcy, and well materials.

| Family | Main fields and component counts |
|---|---|
| H1 Darcy | `Pressure`/`Solution`: 1; `Derivative`/`GradU`: fDim; `Flux`/`MinusKGradU`: 3; `POrder`: 1 |
| Mixed Darcy | `Flux`: 3; `Pressure`: 1; `GradPressure`: 3; `DivFlux`/`Divergence`: 1; `POrder`: 1 |
| H1 well | `Pressure`: 1; `GradPressure`: fDim; `Flux`: fDim (normally 1); `POrder`: 1 |
| Mixed well | `Flux`, `Pressure`, `DivFlux`/`Divergence`, `GradPressure`: 1 each |

For H1 Darcy, `KDuDx/y/z` are components of `K grad(p)` and `Flux` has the opposite
sign. H1 Darcy `Errors()` fills three squared pointwise quantities: pressure plus
permeability-weighted gradient error, pressure L2 error, and the weighted gradient
error. It expects an exact pressure/gradient callback and a pre-sized three-entry
output. NeoPZ performs integration; these pointwise values are not final norms.

Mixed Darcy `Errors()` resizes and fills five squared pointwise quantities:
pressure, flux, divergence residual, pressure gradient, and flux-plus-divergence.
Its exact callback supplies pressure and gradient; it reconstructs exact flux
with `-K grad(p)` and obtains exact divergence from the forcing callback.

Several historical output names are misleading or incomplete:

- H1 Darcy and H1 well `Divergence` currently sum components of a reconstructed
  flux/gradient quantity; they do not differentiate flux to obtain its divergence.
- H1 Darcy and H1 well `ExactFlux` copy the derivative matrix of the exact callback
  directly, while mixed Darcy transforms an exact gradient with `-K`.
- `TPZWannDarcyNL::Solution()` implements only IDs 1, 2, 7, 8, 9, and 10, despite
  inheriting a broader name registry.
- Linear H1 Darcy and H1 well register error/`FluxL2` fields without corresponding
  branches in their current `Solution()` implementations.
- Mixed Darcy's `Permeability` field currently returns numerical flux divergence;
  `normalFlux` returns the first flux component. `ExactDiv` reads derivative entry
  `(2,0)` from its exact callback, while `ExactDivSigma` uses the forcing callback.
- Mixed Darcy's `g_average`/`u_average` outputs require data entries 2/3. Error
  visualization names are registered but have no populated output branches here.

These are descriptions of the existing implementation, not changes to its
numerical behavior. See [well specifics](nonlinear-wells.md) before interpreting
well output or exact-error diagnostics.
