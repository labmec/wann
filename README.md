# WANN paper experiments

This branch contains the three numerical experiments from the wellbore–reservoir coupling paper.

## Build

Use a C++17 compiler, CMake 3.14 or newer, and an installed NeoPZ build. From the repository root, run

```sh
cmake -S . -B build
cmake --build build -j 4
```

If CMake cannot find NeoPZ, configure with `-DNeoPZ_DIR=/path/to/lib/cmake/neopz`. The targets prefer MKL, then MUMPS, then Skyline, according to the installed NeoPZ configuration.

## Run the experiments

Run these commands from the repository root. Each target writes its output in the directory where it runs. Input filenames are resolved relative to `input/`.

```sh
mkdir -p runs/radial runs/adaptive runs/coupled
```

### Section 6.1 — Radial flow

```sh
(cd runs/radial && ../../build/targets/radialFlowTest radialTest.json 4)
```

The target runs the mixed degree-1 and H1 degree-1 and degree-3 approximations on mesh levels 0 through 4. It writes the DoF counts and relative errors to `runs/radial/radialTestResults.txt`, and exports the final solution for each formulation as VTK. The command-line refinement limit defaults to 2. Use `0` for a coarse-mesh check.

### Section 6.2 — Adaptive refinement

```sh
(cd runs/adaptive && ../../build/targets/adaptiveRefinement ozkan1999adaptive.json 3 > adaptive.log)
wc -l runs/adaptive/ozkan1999_refProcess.txt
```

The target solves two reservoir Darcy problems on each mesh, reports total and maximum element errors in scaled and physical units, and writes the mesh and solution VTKs for each iteration. The command-line limit defaults to one pass. With `3`, the target allows three passes and can stop earlier if the scaled estimated error meets the tolerance or the mixed system reaches the equation limit. The refinement history contains one line per completed pass. With post-processing verbosity enabled, `PragerSyngeEst_0.vtk`, `PragerSyngeEst_1.vtk`, and subsequent files contain the elementwise error field for each mesh.

### Section 6.3 — Coupled Ozkan case

For the paper configuration, check that the adaptive run produced three history lines. Copy that history next to the bundled Ozkan mesh, then make a refined copy of the coupled JSON.

```sh
cp runs/adaptive/ozkan1999_refProcess.txt input/ozkan1999_refProcess.txt
cp input/ozkan1999.json input/ozkan1999refined.json
```

Set `MeshData.customRefinement` to `1` and `MeshData.NumUniformRef` to `0` in `input/ozkan1999refined.json`, then run

```sh
(cd runs/coupled && ../../build/targets/wann3dNL ozkan1999refined.json)
```

The bundled `ozkan1999.json` applies one uniform refinement. Custom refinement takes precedence over uniform refinement. The bundled history contains two adaptive passes and can be replayed with the same JSON changes above. A refinement history must come from the same starting mesh and geometry settings.

The target runs the mixed and H1 coupled solves and writes VTK files for both. Check the Newton convergence messages for each solve. The three CSV files used for the mixed well comparison are `fluxOzkanWell.csv` (STB/day), `divergenceOzkanWell.csv` (dimensionless inflow), and `pressureOzkanWell.csv` (psi pressure drop). All three use normalized distance from heel to toe.

The experiments use the meshes bundled in `input/`. See the [Gmsh instructions](meshing/README.md) to regenerate them and the [C++ documentation](docs/README.md) for the implementation.
