Include "params.geo";

// Import both meshes as discrete entities
Merge "reservoir.msh";
Merge "nearWell.msh";

// Merge duplicated nodes that are coincident within tolerance
Coherence Mesh;

// Set physical groups for the combined mesh
// A bit harcoded here. How to do it cleaner?
Physical Volume("volume_reservoir",310) = Volume "*";
Physical Surface("surface_cap_rock",306) = {43,44,11,20};

// Write merged mesh
Save "../input/combined.msh";
