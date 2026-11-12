# WANN paper experiments

This branch contains the three numerical experiments for the wellbore–reservoir coupling paper. It requires a C++17 compiler, CMake, and an installed NeoPZ build.

```sh
cmake -S . -B build
cmake --build build -j 4
```

If CMake cannot find NeoPZ, pass `-DNeoPZ_DIR=/path/to/lib/cmake/neopz` when configuring.

## Run the experiments

Run the executables from a separate output directory. The JSON and mesh files are read from `input/`, and results are written to the current directory.

```sh
mkdir -p runs
cd runs
../build/targets/radialFlowTest radialTest.json 4
../build/targets/adaptiveRefinement ozkan1999adaptive.json 3
```

- **Section 6.1, radial flow.** `radialTestResults.txt` contains the DoF counts and relative errors for the mixed and both H1 approximations. Passing `4` runs levels 0 through 4.
- **Section 6.2, adaptive refinement.** The target writes the error history, VTK files, and `ozkan1999_refProcess.txt`. Passing `3` allows three refinement passes.
- **Section 6.3, coupled Ozkan case.** To use the mesh from Section 6.2, copy the generated `ozkan1999_refProcess.txt` beside `ozkan1999.msh` in `input/` and set `MeshData.customRefinement` to `1` in `input/ozkan1999.json`. Then run `../build/targets/wann3dNL ozkan1999.json` from the output directory. The target writes VTK files and the three mixed well profiles as CSVs.

The checked-in Ozkan refinement history contains two passes. Regenerate it for a three-pass run. A refinement history belongs to the specific mesh that generated it.

The experiments use the meshes bundled in `input/`. See the [Gmsh instructions](meshing/README.md) to regenerate the Ozkan and radial meshes. The [C++ documentation](docs/README.md) covers the implementation.
