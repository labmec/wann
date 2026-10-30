Include "params.geo";

// Import both meshes as discrete entities
Merge "reservoir.msh";
Merge "nearWell.msh";

// Merge duplicated nodes that are coincident within tolerance
Coherence Mesh;
Save "combined.msh";
