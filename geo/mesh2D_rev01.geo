SetFactory("OpenCASCADE");

Lw = 1; //[m]
Dw = 0.1; //[m]

Lr = 10; //[m]
Wr = 5; //[m]
Hr = 2; //[m] 

h_farfield = Hr/10.;
h_wellbore = Dw/5;

Rectangle(1) = {0., -Dw/2, 0, Lw, Dw, 0};
//Rectangle(2) = {-(Lr-Lw)/2, -Wr/2, 0, Lr, Wr, 0};
//all_surf() = BooleanDifference{ Surface{2}; Delete; }{ Surface{1}; };

Point(5) = {-(Lr-Lw)/2, -Wr/2, 0, h_farfield};
Point(6) = {( Lr+Lw)/2, -Wr/2, 0, h_farfield};
Point(7) = {-(Lr-Lw)/2,  Wr/2, 0, h_farfield};
Point(8) = {( Lr+Lw)/2,  Wr/2, 0, h_farfield};

Line(5) = {5,6};
Line(6) = {5,7};
Line(7) = {6,8};
Line(8) = {7,8};

Line(9)  = {5,1};
Line(10) = {6,2};
Line(11) = {8,3};
Line(12) = {7,4};

Curve Loop(2) = {9,4,-12,-6};
Curve Loop(3) = {5,10,1,-9};
Curve Loop(4) = {-10,7,11,-2};
Curve Loop(5) = {12,3,-11,-8};

Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

//Line{9,10,11,12} In Surface{2};

v_div = 20;
h_div = 60;
i_div = 50;
Transfinite Curve{2,4} = v_div Using Bump 0.3;
Transfinite Curve{1,3} = h_div Using Bump 0.3;
Transfinite Curve{6,7} = v_div;
Transfinite Curve{5,8} = h_div;
Transfinite Curve{9,10,11,12} = i_div Using Progression 0.9;
Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Transfinite Surface{4};
Transfinite Surface{5};

//MeshSize {5, 6, 7, 8} = h_farfield; //far field
//MeshSize {1, 2, 3, 4} = h_wellbore; //wellbore

Physical Surface(111) = {1}; //wellbore
Physical Surface(112) = {2,3,4,5}; //reservoir

Physical Line(100) = {5,6,7,8};//reservoir farfield
Physical Line(101) = {4};//wellbore heel
Physical Line(102) = {2};//wellbore toe
Physical Line(103) = {1,3};//wellbore cylinder

General.NumThreads = 4;
Mesh.OptimizeThreshold = 0.5;
Mesh.Smoothing = 10;
//Mesh.Algorithm = 8; //(1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh.Algorithm = 6; //(1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh 2;
RecombineMesh;
