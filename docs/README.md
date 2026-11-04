# WANN C++ documentation

WANN couples stationary reservoir Darcy flow to one-dimensional well flow and
exports results for analysis and neural-network training. This documentation
covers the project-owned C++ library in `src/`, its simulation drivers, and the
historical adaptivity code. The bundled `src/json.hpp` and NeoPZ classes are
external dependencies rather than WANN APIs.

Start with [architecture](architecture.md) to understand the responsibilities
and data flow. Then read the relevant implementation guide:

- [Configuration, meshes, and solver workflow](workflow.md)
- [Darcy materials, tensor dimensions, and outputs](materials.md)
- [Well friction and Newton assembly](nonlinear-wells.md)
- [Export and well diagnostics](postprocessing.md)
- [Legacy code and current limitations](legacy.md)
- [Reservoir transport: first increment](transport.md)

The headers contain Doxygen comments for each project-owned class and method,
including private mesh-building helpers. Implementation comments explain selected
mathematical and storage conventions. These documents describe the current code;
API declarations or input options do not imply that every historical path works.

## Optional browsable API reference

Doxygen is optional and is not part of the simulation build. From the repository
root, with Doxygen already installed, run:

```sh
doxygen docs/Doxyfile
```

The HTML entry point is `docs/api/html/index.html`. Generated files are local
build products and should not be committed. The configuration includes private
helpers and the historical adaptivity header, excludes bundled JSON and generated
path headers, and does not require Graphviz or LaTeX. It does not modify the CMake
build or enable historical source files.

## Maintaining the documentation

Keep a method's contract next to its declaration. Describe its purpose, required
input state, changed outputs, units/frame conventions, and exceptional or
unimplemented behavior. Use `@param[in,out]` for additive assembly matrices and
vectors: `Contribute()` does not clear them. Describe relationships between classes
in these guides rather than duplicating long explanations in every header.

When changing an implementation later, update the corresponding contracts and
limitations at the same time. Distinguish an intended mathematical quantity from
the quantity the code currently computes, especially for post-processing fields
and experimental residual/line-search routines.
