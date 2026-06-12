Include "params.geo";

s = Sin(Pi/4.);

// =============================================================================

// Centers for circular arcs
c0 = newp; Point(c0) = {Lw/2, 0, -Hr/2, 1.0};

// External points
p1 = newp; Point(p1) = {Lw/2 + Wr, 0, -Hr/2, 1.0};
p2 = newp; Point(p2) = {Lw/2, Wr, -Hr/2, 1.0};
p3 = newp; Point(p3) = {Lw/2 - Wr, 0, -Hr/2, 1.0};
p4 = newp; Point(p4) = {Lw/2, -Wr, -Hr/2, 1.0};

// Internal points
p5 = newp; Point(p5) = {-Hr/2, -Hr/2, -Hr/2, 1.0};
p6 = newp; Point(p6) = {Lw+Hr/2, -Hr/2, -Hr/2, 1.0};
p7 = newp; Point(p7) = {Lw+Hr/2, Hr/2, -Hr/2, 1.0};
p8 = newp; Point(p8) = {-Hr/2, Hr/2, -Hr/2, 1.0};

// Create outer ball
l1 = newl; Circle(l1) = {p1, c0, p2};
l2 = newl; Circle(l2) = {p2, c0, p3};
l3 = newl; Circle(l3) = {p3, c0, p4};
l4 = newl; Circle(l4) = {p4, c0, p1};

// Create inner square (SmallBox)
l5 = newl; Line(l5) = {p5, p6};
l6 = newl; Line(l6) = {p6, p7};
l7 = newl; Line(l7) = {p7, p8};
l8 = newl; Line(l8) = {p8, p5};

Transfinite Line {l5, l7} = axial_div;
Transfinite Line {l6, l8} = h_div+1;

// Grade mesh size with distance from the inner square (SmallBox):
// keep fine elements near l5..l8 and grow toward the outer boundary.
Field[1] = Distance;
Field[1].CurvesList = {l5, l6, l7, l8};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = size_min_res;
Field[2].SizeMax = size_max_res;
Field[2].DistMin = dist_min_res;
Field[2].DistMax = dist_max_res;

Background Field = 2;

// Use only the background field for size control in the domain interior.
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

// Build a planar surface with a hole: BigBox minus SmallBox
clOuter = newcl; Curve Loop(clOuter) = {l1, l2, l3, l4};
clInner = newcl; Curve Loop(clInner) = {l5, l6, l7, l8};
sf = news; Plane Surface(sf) = {clOuter, clInner};

// Extrude the resulting surface in z direction
v[] = Extrude {0, 0, Hr} {Surface{sf}; Layers{h_div}; Recombine;};

// Physical groups for the final mesh
Physical Surface("surface_farfield",307) = {12, 13, 14, 15};

// Physical groups for merging with the near-well mesh
Physical Surface("surface_cap_rock",306) = {11, 20};
Physical Volume("volume_reservoir",310) = {1};

// Generate mesh
Mesh 3;
Save "reservoir.msh";
