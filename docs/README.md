# WANN C++ documentation

WANN couples stationary reservoir Darcy flow to one-dimensional well flow and
exports results for the paper experiments. This documentation
covers the project-owned C++ library in `src/`, its simulation drivers, and the
adaptive-refinement code. The bundled `src/json.hpp` and NeoPZ classes are
external dependencies rather than WANN APIs.

For build and run commands, see the [paper experiment instructions](../README.md).
Start with [architecture](architecture.md) for the code structure and data flow. Then read the relevant implementation guide:

- [Configuration, meshes, and solver workflow](workflow.md)
- [Darcy materials, tensor dimensions, and outputs](materials.md)
- [Well friction and Newton assembly](nonlinear-wells.md)
- [Export and well diagnostics](postprocessing.md)
- [Paper targets and implementation limits](legacy.md)

The headers contain Doxygen comments for each project-owned class and method,
including private mesh-building helpers. Implementation comments explain selected
mathematical and storage conventions.

## Optional browsable API reference

Doxygen is optional and is not part of the simulation build. From the repository
root, with Doxygen already installed, run:

```sh
doxygen docs/Doxyfile
```

The HTML entry point is `docs/api/html/index.html`. Generated files are local
build products and should not be committed. The configuration includes private
helpers and the adaptivity header, excludes bundled JSON and generated
path headers, and does not require Graphviz or LaTeX.

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
