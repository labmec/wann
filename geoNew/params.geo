SetFactory("OpenCASCADE"); // Geometry kernel

// Shared parameters used by both geometry generators.
Geometry.MatchMeshTolerance = 1e-12;
Geometry.Tolerance = 1e-16;
Mesh.ToleranceReferenceElement = 1e-10;
General.Terminal = 1;
Geometry.SnapPoints = 0;

// Well bore dimensions
Lw = 200;
Rw = 0.1;
ecc = -0.5; // Eccentricity of the wellbore (-0.5 <= ecc <= 0.5)

// Reservoir dimensions
Hr = 20;
Wr = 300;
Lr = 400;

// Meshing parameters shared by both geometries
h_div = 5;        // Division in z
axial_div = 16;   // Divisions along the well length

// Well only meshing parameters
radial_div = 8;   // Radial division of the near-well region
p_res = 1.6;      // Progression coefficient of the near-well region mesh
p_well = 0.3;     // Progression coefficient of the wellbore mesh

// Reservoir only meshing parameters
size_min_res = 8.0;
size_max_res = 80.0;
dist_min_res = 10.0;
dist_max_res = 100.0;
