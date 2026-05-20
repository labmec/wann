SetFactory("OpenCASCADE"); // Geometry kernel

// Tolerance parameters
Geometry.MatchMeshTolerance = 1e-12;
Geometry.Tolerance = 1e-16;
Mesh.ToleranceReferenceElement = 1e-10;
General.Terminal = 1;
Geometry.SnapPoints = 0;

// Reservoir dimensions
Lr = 400;
Hr = 20;
Wr = 300;

// Well length
Lw = 200;

h_div = 5; // Division in z
axial_div = 16; // Divisions along the well length

res_length_div = 6; // Divisions along the reservoir length
res_width_div = 4; // Divisions along the reservoir width
diagonal_div = 10; // Divisions along the diagonal lines
prog_diag = 0.8; // Progression coefficient of the diagonal lines

s = Sin(Pi/4.);

// =============================================================================

// External points
p1 = newp; Point(p1) = {-(Lr-Lw)/2, -Wr/2, -Hr/2, 1.0};
p2 = newp; Point(p2) = {(Lr+Lw)/2, -Wr/2, -Hr/2, 1.0};
p3 = newp; Point(p3) = {(Lr+Lw)/2, Wr/2, -Hr/2, 1.0};
p4 = newp; Point(p4) = {-(Lr-Lw)/2, Wr/2, -Hr/2, 1.0};

// Internal points
p5 = newp; Point(p5) = {-Hr/2, -Hr/2, -Hr/2, 1.0};
p6 = newp; Point(p6) = {Lw+Hr/2, -Hr/2, -Hr/2, 1.0};
p7 = newp; Point(p7) = {Lw+Hr/2, Hr/2, -Hr/2, 1.0};
p8 = newp; Point(p8) = {-Hr/2, Hr/2, -Hr/2, 1.0};

// Create outer square (BigBox)
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

// Create inner square (SmallBox)
l5 = newl; Line(l5) = {p5, p6};
l6 = newl; Line(l6) = {p6, p7};
l7 = newl; Line(l7) = {p7, p8};
l8 = newl; Line(l8) = {p8, p5};

// Diagonal lines
l9 = newl; Line(l9) = {p1, p5};
l10 = newl; Line(l10) = {p2, p6};
l11 = newl; Line(l11) = {p3, p7};
l12 = newl; Line(l12) = {p4, p8};

Transfinite Line {l5, l7} = axial_div;
Transfinite Line {l6, l8} = h_div;
Transfinite Line {l1, l3} = res_length_div;
Transfinite Line {l2, l4} = res_width_div;
Transfinite Line {l9, l10, l11, l12} = diagonal_div Using Progression prog_diag;

// Surfaces
cl1 = newcl; Curve Loop(cl1) = {l9, -l8, -l12, l4};
cl2 = newcl; Curve Loop(cl2) = {l1, l10, -l5, -l9};
cl3 = newcl; Curve Loop(cl3) = {l2, l11, -l6, -l10};
cl4 = newcl; Curve Loop(cl4) = {-l11, l3, l12, -l7};

sf1 = news; Plane Surface(sf1) = {cl1};
sf2 = news; Plane Surface(sf2) = {cl2};
sf3 = news; Plane Surface(sf3) = {cl3};
sf4 = news; Plane Surface(sf4) = {cl4};

// Extrude surfaces in z direction
v1[] = Extrude {0, 0, Hr} {Surface{sf1}; Layers{h_div}; Recombine;};
v2[] = Extrude {0, 0, Hr} {Surface{sf2}; Layers{h_div}; Recombine;};
v3[] = Extrude {0, 0, Hr} {Surface{sf3}; Layers{h_div}; Recombine;};
v4[] = Extrude {0, 0, Hr} {Surface{sf4}; Layers{h_div}; Recombine;};

// Generate mesh
Mesh 3;
Save "meshTestA.msh";

