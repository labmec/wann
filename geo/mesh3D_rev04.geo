SetFactory("OpenCASCADE");

//wellbore dimensions
Lw = 1; //[m]
Dw = 0.1; //[m]

//block dimensions
Lb = 2*Lw; //[m]
Wb = 8*Dw; //[m]
Hb = 8*Dw; //[m] 

//reservoir dimensions
Lr = 10; //[m]
Wr = 8; //[m]
Hr = 6; //[m]

//wellbore division
h_div = 6;//horizontal division
v_div = 3;//vertical division
w_div = 3;//width division

//reservoir division (only one parameter for now)
h_farfield = Hr/10.;

Box(1) = {0, -Dw/2, -Dw/2, Lw, Dw, Dw};//wellbore
Box(2) = {-(Lb-Lw)/2, -Wb/2, -Hb/2, Lb, Wb, Hb};//block
//remove wellbore volume
all_vol() = BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; };
//all_vol() = BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }; //used to generate the block+wellbore mesh only

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

Box(9) = {-(Lr-Lw)/2, -Wr/2, -Hr/2, Lr, Wr, Hr};//block
all_vol() = BooleanDifference{ Volume{3,4,5,6,7,8,9}; }{ Volume{1,2}; Delete; };
Delete{ Volume{9}; }
//Delete{ Volume{2};} //used to generate the block+wellbore mesh only

//Transfinite Curve{:} = 2;
Transfinite Curve{9,10,11,12,17,19,20,22} = h_div Using Bump 0.3;
Transfinite Curve{25,26,27,28,29,30,31,32} = v_div Using Progression 1.2;
Transfinite Curve{1,2,3,4,5,6,7,8,13,14,15,16,18,21,23,24} = w_div Using Progression 1;
//Transfinite Surface{:};
Transfinite Surface{1,2,3,4,5,6,7,8,9,10,11,12,100,101,102,103,104,105,106,107,108,109,110,111};
Transfinite Volume{3,4,5,6,7,8};
//Recombine Surface{:};
//Recombine Surface{1,2,3,4,5,6,7,8,9,10,11,12,100,101,102,103,104,105,106,107,108,109,110,111};

//far field nodes
//MeshSize {17, 18, 19, 20, 21, 22, 23, 24} = h_farfield; //far field

Physical Curve("curve_wellbore",100) = {9,10,11,12};
Physical Curve("curve_toe",101) = {5,6,7,8};
Physical Curve("curve_heel",102) = {1,2,3,4};
Physical Surface("surface_wellbore_cylinder",103) = {3,4,5,6};//wellbore cylinder
Physical Surface("surface_wellbore_toe",104) = {2};//wellbore toe
Physical Surface("surface_wellbore_heel",105) = {1};//wellbore heel
Physical Surface("surface_farfield",106) = {112,113,114,115};//reservoir farfield
Physical Volume("volume_reservoir",107) = {10}; //reservoir

//Mesh.SubdivisionAlgorithm=2;  //all hexas
General.NumThreads = 4;
Mesh.OptimizeThreshold = 0.5;
Mesh.Smoothing = 40;
Mesh.Algorithm = 8;//1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh.Algorithm3D = 1; //1: Delaunay (default), 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT (paralelizável)
//Mesh.SubdivisionAlgorithm = 2; //(0: none, 1: all quadrangles, 2: all hexahedra, 3: barycentric)
//Mesh.RecombinationAlgorithm = 3; //(0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)
Mesh 3;
OptimizeMesh "Netgen";
OptimizeMesh "Gmsh"; // "Netgen" is slow

