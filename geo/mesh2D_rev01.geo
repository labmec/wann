SetFactory("OpenCASCADE");

Lw = 1; //[m]
Dw = 0.1; //[m]

Lr = 10; //[m]
Wr = 5; //[m]
Hr = 2; //[m] 

h_farfield = Hr/10.;
h_wellbore = Dw/5;

Rectangle(1) = {0., -Dw/2, 0, Lw, Dw, 0};
Rectangle(2) = {-(Lr-Lw)/2, -Wr/2, 0, Lr, Wr, 0};
//Box(2) = {-(Lr-Lw)/2, -Wr/2, -Hr/2, Lr, Wr, Hr};

all_surf() = BooleanDifference{ Surface{2}; Delete; }{ Surface{1}; };

MeshSize {5, 6, 7, 8} = h_farfield; //far field
MeshSize {1, 2, 3, 4} = h_wellbore; //wellbore

Physical Surface(111) = {1}; //wellbore
Physical Surface(112) = {2}; //reservoir

Physical Line(100) = {5,6,7,8};//reservoir farfield
Physical Line(101) = {4};//wellbore heel
Physical Line(102) = {2};//wellbore toe
Physical Line(103) = {1,3};//wellbore cylinder

General.NumThreads = 4;
Mesh.OptimizeThreshold = 0.5;
Mesh.Smoothing = 10;
Mesh.Algorithm = 8; //(1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh 2;
RecombineMesh;