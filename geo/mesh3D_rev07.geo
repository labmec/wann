SetFactory("OpenCASCADE");

Geometry.MatchMeshTolerance = 1e-10;
Geometry.Tolerance = 1e-12;
Mesh.ToleranceReferenceElement = 1e-10;

//wellbore dimensions
Lw = 1; //[m] wellbore length
Dw = 0.1; //[m] wellbore diameter

//reservoir dimensions
Lr = 4; //[m] length
Wr = 3; //[m] width
Hr = 1; //[m] height

//wellbore division
h_div = 8;  //horizontal division of the reservoir
r_div = 4;  //radius division (of each circle quarter)
l_div = 10; //axial division of the wellbore
p_res = 1.6;//progression coeficient of the reservoir mesh
p_well = 0.3;//progression coeficient of the wellbore mesh (bump scheme)

s = Sin(Pi/4.);

//points of the heel part {{{
Point(1) = {0, 0, 0, 1.0};//center of the domain (wellbore heel)
Point(2) = {0, -Dw*s, Dw*s, 1.0};
Point(3) = {0, Dw*s, Dw*s, 1.0};
Point(4) = {0, Dw*s, -Dw*s, 1.0};
Point(5) = {0, -Dw*s, -Dw*s, 1.0};
Point(6)= {-(Lr-Lw)/2, -Wr/2, Hr/2, 1.0};//FIXME
Point(7)= {-(Lr-Lw)/2, Wr/2, Hr/2, 1.0};
Point(8)= {-(Lr-Lw)/2, Wr/2, -Hr/2, 1.0};
Point(9)= {-(Lr-Lw)/2, -Wr/2, -Hr/2, 1.0};
//}}}
//points of the toe part {{{
Point(10) = {Lw, 0, 0, 1.0};//center of the domain (wellbore toe)
Point(11) = {Lw, -Dw*s, Dw*s, 1.0};
Point(12) = {Lw, Dw*s, Dw*s, 1.0};
Point(13) = {Lw, Dw*s, -Dw*s, 1.0};
Point(14) = {Lw, -Dw*s, -Dw*s, 1.0};
Point(15)= {Lw+(Lr-Lw)/2, -Wr/2, Hr/2, 1.0};
Point(16)= {Lw+(Lr-Lw)/2, Wr/2, Hr/2, 1.0};
Point(17)= {Lw+(Lr-Lw)/2, Wr/2, -Hr/2, 1.0};
Point(18)= {Lw+(Lr-Lw)/2, -Wr/2, -Hr/2, 1.0};
//}}}
//circles and lines of the wellbore (heel part) {{{
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line(5) = {7, 3};
Line(6) = {8, 4};
Line(7) = {9, 5};
Line(8) = {6, 2};
Line(9) = {7, 6};
Line(10) = {6, 9};
Line(11) = {9, 8};
Line(12) = {8, 7};
//}}}
//circles and lines of the wellbore (toe part) {{{
Circle(13) = {11, 10, 12};
Circle(14) = {12, 10, 13};
Circle(15) = {13, 10, 14};
Circle(16) = {14, 10, 11};
Line(17) = {11, 15};
Line(18) = {12, 16};
Line(19) = {13, 17};
Line(20) = {14, 18};
Line(21) = {15, 16};
Line(22) = {16, 17};
Line(23) = {17, 18};
Line(24) = {18, 15};
//}}}
//curve loops and surfaces of the heel part {{{
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, -1, -8, -9};
Surface(2) = {2};
Curve Loop(4) = {10, 7, 4, -8};
Surface(3) = {4};
Curve Loop(6) = {7, -3, -6, -11};
Surface(4) = {6};
Curve Loop(8) = {6, -2, -5, -12};
Surface(5) = {8};
Curve Loop(10) = {12, 9, 10, 11};
Plane Surface(6) = {10};
Surface Loop(1) = {2, 5, 4, 6, 3, 1};
Volume(1) = {1};
//}}}
//curve loops and surfaces of the heel part {{{
Curve Loop(12) = {13, 14, 15, 16};
Surface(7) = {12};
Curve Loop(14) = {13, 18, -21, -17};
Surface(8) = {14};
Curve Loop(16) = {14, 19, -22, -18};
Surface(9) = {16};
Curve Loop(18) = {15, 20, -23, -19};
Surface(10) = {18};
Curve Loop(20) = {16, 17, -24, -20};
Surface(11) = {20};
Curve Loop(22) = {24, 21, 22, 23};
Plane Surface(12) = {22};
Surface Loop(2) = {8, 7, 11, 10, 9, 12};
Volume(2) = {2};
//}}}
//other lines and surfaces {{{
Line(25) = {2, 11};
Line(26) = {3, 12};
Line(27) = {4, 13};
Line(28) = {5, 14};
Line(29) = {6, 15};
Line(30) = {7, 16};
Line(31) = {8, 17};
Line(32) = {9, 18};
Curve Loop(24) = {9, 29, 21, -30};
Plane Surface(13) = {24};
Curve Loop(26) = {30, 22, -31, 12};
Plane Surface(14) = {26};
Curve Loop(28) = {31, 23, -32, 11};
Plane Surface(15) = {28};
Curve Loop(30) = {32, 24, -29, 10};
Plane Surface(16) = {30};

Curve Loop(32) = {4, 25, -16, -28};
Surface(17) = {32};
Curve Loop(34) = {1, 26, -13, -25};
Surface(18) = {34};
Curve Loop(36) = {2, 27, -14, -26};
Surface(19) = {36};
Curve Loop(38) = {27, 15, -28, -3};
Surface(20) = {38};
Curve Loop(40) = {5, 26, 18, -30};
Plane Surface(21) = {40};
Curve Loop(42) = {8, 25, 17, -29};
Plane Surface(22) = {42};
Curve Loop(44) = {27, 19, -31, 6};
Plane Surface(23) = {44};
Curve Loop(46) = {28, 20, -32, 7};
Plane Surface(24) = {46};

Surface Loop(3) = {13, 2, 21, 8, 22, 18};
Volume(3) = {3};
Surface Loop(4) = {21, 14, 5, 9, 23, 19};
Volume(4) = {4};
Surface Loop(5) = {15, 23, 4, 10, 24, 20};
Volume(5) = {5};
Surface Loop(6) = {16, 3, 11, 24, 22, 17};
Volume(6) = {6};
//}}}
//apply transfinite progression to lines {{{
Transfinite Line {-5, -6, -7, -8} = h_div Using Progression p_res; //reservoir
Transfinite Line {17, 18, 19, 20} = h_div Using Progression p_res; //reservoir
Transfinite Line {1, 2, 3, 4, 9, 10, 11, 12} = r_div; //wellbore radius 
Transfinite Line {13, 14, 15, 16, 21, 22, 23, 24} = r_div; //wellbore radius
Transfinite Line {25, 26, 27, 28} = l_div Using Bump p_well; //wellbore length (with bump)
Transfinite Line {29, 30, 31, 32} = l_div; //horizotanl reservoir lines (no bump)

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
///}}}
//set physical entites {{{
Physical Curve("curve_wellbore",100) = {26};//just one line of the wellbore
Physical Curve("curve_toe",101) = {13,14,15,16};
Physical Curve("curve_heel",102) = {1,2,3,4};
Physical Surface("surface_wellbore_cylinder",103) = {17,18,19,20};//wellbore cylinder
Physical Surface("surface_wellbore_toe",104) = {7};//wellbore toe
Physical Surface("surface_wellbore_heel",105) = {1};//wellbore heel
Physical Surface("surface_farfield",106) = {6,12,14,16};//reservoir farfield
Physical Surface("surface_cap_rock",107) = {13,15};//cap rock
Physical Point("point_heel",108) = {3};
Physical Point("point_toe",109) = {12};
all_volumes[] = Volume "*";
Physical Volume("volume_reservoir",108) = all_volumes[]; //reservoir
//}}}

Delete{ Point{1}; Point{10};
}

//General.NumThreads = 4;
//Mesh.OptimizeThreshold = 0.5;
//Mesh.Smoothing = 40;
Mesh.Smoothing = 0;//smoothing is changing the coordinates of the cylinder (which is not desirable here)
//Mesh.Algorithm = 8;
//Mesh 3;
//OptimizeMesh "Netgen";
//OptimizeMesh "Gmsh"; // "Netgen" is slow


