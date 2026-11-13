# Well friction and Newton assembly

The well materials describe one-dimensional axial flow. Their physical constructor
accepts material ID, diameter `Dw`, viscosity `mu`, density `rho`, reservoir
reference pressure `pres`, and pseudo-resistivity `Kvw`. Positive physical
coefficients are required by divisions and fractional powers in the formulas.
A default or dimension-only constructor does not initialize all friction data.

Let

```text
CL = pi Dw^4 / (128 mu)
c  = 2.252610888 Dw^(19/7) / (mu^(1/7) rho^(3/7))
```

The H1 well stores `fCLin = CL` and `fC = c`. The mixed well stores their resistance
counterparts `fCLin = 1/CL` and `fC = c^(-7/4)`.

## Pressure-only H1 well

`TPZWannH1Well` uses the pressure solution and its first derivative. In the laminar
regime the reconstructed flow is `Q = -CL * dp/ds`; turbulent assembly and flux
output use `Q = -sign(dp/ds) * c * abs(dp/ds)^(4/7)`, where `s` is the local
axial coordinate.

Full `Contribute()` computes a critical flow from Reynolds number 1187.38 and
converts it to a critical gradient using `CL`. The first Newton iteration remains
laminar. Later iterations switch when the absolute pressure derivative exceeds
that threshold. The currently installed volume forcing callback is not evaluated in these well
assembly paths. It accumulates tangent and negative-residual contributions,
including `Kvw*(pres-p)` as the optional pseudo-resistivity term.

The turbulent negative residual retains the gradient sign. Its tangent coefficient
is `c * (4/7) * abs(dp/ds)^(-3/7)`, positive for either gradient sign, consistent
with the derivative of the signed flow law.

When `fAssembleRHSOnly` is true, `ContributeResidual()` estimates Reynolds from a
laminar flow estimate, then rechecks the turbulent estimate and stops if it falls
back below the transition threshold. It adds no tangent matrix.

Boundary type 0 contributes nothing because analysis imposes pressure values.
Types 1 and 2 add flow-load/Robin terms. `Solution()` supplies only scalar
`Pressure` and signed axial `Flux`. Flux uses the local axial pressure derivative,
with its sign following the element orientation.

`Errors()` uses unit gradient weighting against an exact callback, not the
nonlinear friction energy. It expects three output entries. The stored error
names are H1, L2, and H1 seminorm.

## Mixed flow/pressure well

`TPZWannMixedWell` uses flux-space entry zero and pressure-space entry one. It
extracts only row zero of `fDeformedDirections` and component zero of the numerical
flow. This is an x-component assumption even when geometry/pressure grouping uses
a more general well axis.

The constitutive resistance assembled by the full path is:

```text
laminar:   (1/CL) Q
turbulent: c^(-7/4) Q abs(Q)^(3/4)
```

The turbulent tangent coefficient is `(7/4) c^(-7/4) abs(Q)^(3/4)`. The full path
uses velocity `4Q/(pi Dw^2)` to compute Reynolds, with transition at 1187.38 and a
laminar first iteration. Divergence-pressure coupling and optional
pseudo-resistivity populate the other blocks.

`ContributeResidual()` is a separate implementation. Although it computes a
Reynolds estimate, it then explicitly assigns `turbulent = false`. Consequently,
its residual is always laminar and can differ from the full assembly residual.
This matters when line searches or Taylor diagnostics select residual-only
assembly.

`SetForcingFunctionDual()` stores a callback added against flux-test divergence;
its `order` argument is currently unused. Both full and residual-only paths include
this optional contribution.

Pressure endpoint BCs add a load to flux rows. Prescribed flow endpoint equations
are initialized/filtered by analysis, and type 1 adds no penalty here. Outputs
`Flux`, `Pressure`, `DivFlux`/`Divergence`, and `GradPressure` each have one component.
`Errors()` is unimplemented and calls `DebugStop()`, even though `NEvalErrors()`
returns five.

## Shared flags and Newton diagnostics

Each well class has independent static flags:

- `fIsFirstIteration` forces the initial laminar regime for its full assembly.
- `fAssembleRHSOnly` routes assembly to its residual-only implementation.

The Newton driver resets both first-iteration flags at the start and clears them
after its first correction. Flags are shared across all material instances of a
class, including different wells or analyses.

`LineSearchStep()` supports a full step, bisection, and golden-ratio search.
`ESECANT` is declared but stops as unimplemented. The first Newton correction
always takes a full step regardless of the selected method.

Bisection searches for a root of `dsol^T RHS` and assumes a multiphysics mesh.
Golden-ratio search attempts to reduce the RHS norm over correction scales in
`[0,1]`. Both load trial solutions, reassemble, and toggle only the **mixed-well**
residual-only flag. Bisection's early return leaves that flag enabled. These
routines do not restore the original mesh/RHS state; the Newton loop subsequently
loads its accepted update.

`TaylorTest()` estimates Jacobian consistency by comparing residuals at two
perturbation sizes. It requires an assembled Jacobian/RHS and no line search;
the current implementation is not compatible with equation filtering. It also
uses the mixed-well residual-only flag and leaves a trial solution loaded.
`PostProcessIteration()` provides per-iteration VTK output for debugging, but its
call in the Newton loop is currently commented out.

## Copy behavior

The existing copy/assignment implementations do not fully copy well coefficients.
H1 well assignment copies the base and dimension but omits the physical/friction
members. Mixed well assignment copies most physical members and the turbulent
coefficient, but omits `fCLin` and the dual callback. `NewMaterial()` invokes the
copy constructor, so cloning should not be assumed to preserve all well state.
This documentation records that limitation without changing the copy logic.
