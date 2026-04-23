SetFactory("OpenCASCADE"); // Geometry kernel

// Tolerance parameters
Geometry.MatchMeshTolerance = 1e-12;
Geometry.Tolerance = 1e-16;
Mesh.ToleranceReferenceElement = 1e-10;
General.Terminal = 1;
Geometry.SnapPoints = 0;

SetFactory("OpenCASCADE");

// ----------------------
// Geometry parameters
// ----------------------
Rw = 0.2;
Rr = 1.0;
Lw = 2.0;

Nr = 5;
Nt = 20;   // multiple of 4
Nx = 10;

RadialProgression = 1.3;

// ----------------------
// Geometry definition
// ----------------------

// Inner circle (yz-plane)
Point(1) = {0,  Rw, 0};
Point(2) = {0,  0,  Rw};
Point(3) = {0, -Rw, 0};
Point(4) = {0,  0, -Rw};

// Outer circle
Point(5) = {0,  Rr, 0};
Point(6) = {0,  0,  Rr};
Point(7) = {0, -Rr, 0};
Point(8) = {0,  0, -Rr};

Point(9) = {0,0,0};

// Arcs
Circle(1) = {1,9,2};
Circle(2) = {2,9,3};
Circle(3) = {3,9,4};
Circle(4) = {4,9,1};

Circle(5) = {5,9,6};
Circle(6) = {6,9,7};
Circle(7) = {7,9,8};
Circle(8) = {8,9,5};

// Radial lines
Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

// Surfaces
Curve Loop(1) = {1,10,-5,-9};
Plane Surface(1) = {1};

Curve Loop(2) = {2,11,-6,-10};
Plane Surface(2) = {2};

Curve Loop(3) = {3,12,-7,-11};
Plane Surface(3) = {3};

Curve Loop(4) = {4,9,-8,-12};
Plane Surface(4) = {4};

// ----------------------
// Transfinite definitions
// ----------------------

Transfinite Curve{1,2,3,4,5,6,7,8} = Nt/4;
Transfinite Curve{9,10,11,12} = Nr Using Progression RadialProgression;

Transfinite Surface{1,2,3,4};
Recombine Surface{1,2,3,4};

// ----------------------
// Extrude and CAPTURE volume
// ----------------------
out[] = Extrude {Lw,0,0}
{
  Surface{1,2,3,4};
  Layers{Nx};
  Recombine;
};

Transfinite Volume{out[1], out[3], out[5], out[7]};
Recombine Volume{out[1], out[3], out[5], out[7]};

// ----------------------
// Physical groups
// ----------------------
Physical Volume("volume_reservoir",110) = {1,2,3,4};
Physical Surface("surface_farfield",106) = {19,16,12,7};
Physical Surface("surface_wellbore_cylinder",103) = {18,14,10,5};
Physical Surface("surface_cap_rock",107) = {1,2,3,4,9,13,17,20};
Physical Curve("curve_wellbore",100) = {13};
Physical Point("point_heel",108) = {1};
Physical Point("point_toe",109) = {10};   

// ----------------------
// Meshing
// ----------------------

Delete{ Point{9}; } // Center point no longer needed
Mesh 3;
Coherence Mesh;
Save "../input/radialTest.msh";



