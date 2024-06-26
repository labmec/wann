SetFactory("OpenCASCADE");

Lw = 1; //[m]
Dw = 0.1; //[m]

Lr = 10; //[m]
Wr = 5; //[m]
Hr = 2; //[m] 

h_farfield = Hr/10.;
h_wellbore = Dw/5;

Box(1) = {0, -Dw/2, -Dw/2, Lw, Dw, Dw};
Box(2) = {-(Lr-Lw)/2, -Wr/2, -Hr/2, Lr, Wr, Hr};

all_vol() = BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; };

MeshSize {9, 10, 11, 12, 13, 14, 15, 16}  = h_farfield; //far field
MeshSize {1, 2, 3, 4, 5, 6, 7,8 } = h_wellbore; //wellbore

Physical Volume(111) = {2};
Physical Surface(100) = {7,8,10,12};//reservoir farfield (top and bottom not included)
Physical Surface(101) = {1};//wellbore heel
Physical Surface(102) = {2};//wellbore toe
Physical Surface(103) = {3, 4, 5, 6};//wellbore cylinder

//Mesh.SubdivisionAlgorithm=2;  //all hexas
General.NumThreads = 4;
Mesh.OptimizeThreshold = 0.5;
Mesh.Smoothing = 10;
Mesh.Algorithm3D = 10; //1: Delaunay (default), 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT (paralelizável)

