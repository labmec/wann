SetFactory("OpenCASCADE");

//wellbore dimensions
Lw = 1; //[m]
Dw = 0.1; //[m]

//block dimensions
Lb = 2*Lw; //[m]
Wb = 8*Dw; //[m]
Hb = 8*Dw; //[m] 

//wellbore division (along its axis)
h_div = 10;
v_div = 10;

//h_farfield = Hr/10.;
//h_wellbore = Dw/5;

Box(1) = {0, -Dw/2, -Dw/2, Lw, Dw, Dw};//wellbore
Box(2) = {-(Lb-Lw)/2, -Wb/2, -Hb/2, Lb, Wb, Hb};//block
//remove wellbore volume
all_vol() = BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; };

Line(25) = {1,9};
Line(26) = {2,10};
Line(27) = {3,11};
Line(28) = {4,12};
Line(29) = {6,13};
Line(30) = {5,14};
Line(31) = {7,15};
Line(32) = {8,16};

Curve Loop(100) = {1,25,13,-26};
Curve Loop(101) = {2,27,-14,-25};
Curve Loop(102) = {-3,27,15,-28};
Curve Loop(103) = {4,28,-16,-26};

Curve Loop(104) = {-7,31,-23,-32};
Curve Loop(105) = {-6,30,21,-31};
Curve Loop(106) = {-5,30,-18,-29};
Curve Loop(107) = {-8, 29,24,-32};

Curve Loop(108) = {12,31,-20,-27};
Curve Loop(109) = {10,30,-19,-25};
Curve Loop(110) = {11,32,-22,-28};
Curve Loop(111) = {-9,26,17,-29}; 

Plane Surface(100) = {100};
Plane Surface(101) = {101};
Plane Surface(102) = {102};
Plane Surface(103) = {103};

Plane Surface(104) = {104};
Plane Surface(105) = {105};
Plane Surface(106) = {106};
Plane Surface(107) = {107};

Plane Surface(108) = {108};
Plane Surface(109) = {109};
Plane Surface(110) = {110};
Plane Surface(111) = {111};

Surface Loop(3) = {1,100,101,102,103,7};
Surface Loop(4) = {2,104,105,106,107,12};
Surface Loop(5) = {6,101,109,105,108,9};
Surface Loop(6) = {5,103,107,110,111,11};
Surface Loop(7) = {4,102,104,108,110,10};
Surface Loop(8) = {3,106,109,100,111,8};

Volume(3) = {3};
Volume(4) = {4};
Volume(5) = {5};
Volume(6) = {6};
Volume(7) = {7};
Volume(8) = {8};

all_vol() = BooleanDifference{ Volume{3,4,5,6,7,8}; }{ Volume{2}; Delete; };

Transfinite Curve{:} = 2;
Transfinite Curve{9,10,11,12,17,19,20,22} = h_div Using Bump 0.3;
Transfinite Curve{25,26,27,28,29,30,31,32} = v_div Using Progression 1.2;
Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};

//MeshSize {3, 4, 5, 6, 7, 8, 9, 10} = h_farfield; //far field
//MeshSize {1, 2} = h_wellbore; //wellbore

//Physical Volume(111) = {1}; //wellbore
//Physical Volume(112) = {2}; //reservoir

//Physical Surface(100) = {4,5,7,9};//reservoir farfield
//Physical Surface(101) = {3};//wellbore heel
//Physical Surface(102) = {2};//wellbore toe
//Physical Surface(103) = {1};//wellbore cylinder

//Mesh.SubdivisionAlgorithm=2;  //all hexas
General.NumThreads = 4;
Mesh.OptimizeThreshold = 0.5;
Mesh.Smoothing = 10;
Mesh.Algorithm = 8;//1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh.Algorithm3D = 1; //1: Delaunay (default), 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT (paralelizável)
//Mesh.SubdivisionAlgorithm = 2; //(0: none, 1: all quadrangles, 2: all hexahedra, 3: barycentric)
//Mesh.RecombinationAlgorithm = 3; //(0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)

