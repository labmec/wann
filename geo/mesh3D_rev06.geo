SetFactory("OpenCASCADE");

//wellbore dimensions
Lw = 1; //[m] wellbore length
Dw = 0.1; //[m] wellbore diameter

//reservoir dimensions
Lr = 3; //[m]
Wr = 2; //[m]
Hr = 1; //[m]

//wellbore division
h_div = 6;
t_div = 5; //theta division
r_div = 4; //radius division

//reservoir division (only one parameter for now)
h_farfield = Hr/10.;
i_div = 15; //division of the inclined lines connecting reservoir to wellbore
p_res = 1.2;
b_res = 3.0;

s = Sin(Pi/4.);

Point(1) = {0, 0, 0, 1.0};//center of the domain (wellbore heel)
Point(2) = {0, -0.5*Dw*s, 0.5*Dw*s, 1.0};
Point(3) = {0, -Dw*s, Dw*s, 1.0};
Point(4) = {0, 0.5*Dw*s, 0.5*Dw*s, 1.0};
Point(5) = {0, Dw*s, Dw*s, 1.0};
Point(6) = {0, 0.5*Dw*s, -0.5*Dw*s, 1.0};
Point(7) = {0, Dw*s, -Dw*s, 1.0};
Point(8) = {0, -0.5*Dw*s, -0.5*Dw*s, 1.0};
Point(9) = {0, -Dw*s, -Dw*s, 1.0};
Point(10)= {0, -Wr/2, Hr/2, 1.0}; 
Point(11)= {0, Wr/2, Hr/2, 1.0}; 
Point(12)= {0, Wr/2, -Hr/2, 1.0}; 
Point(13)= {0, -Wr/2, -Hr/2, 1.0}; 

// core lines
Line(1) = {6, 8};
Line(2) = {8, 2};
Line(3) = {2, 4};
Line(4) = {4, 6};
//radial lines of the wellbore
Line(5) = {6, 7};
Line(6) = {8, 9};
Line(7) = {2, 3};
Line(8) = {4, 5};
//circles of the wellbore
Circle(9) = {7, 1, 9};
Circle(10) = {9, 1, 3};
Circle(11) = {3, 1, 5};
Circle(12) = {5, 1, 7};
//reservoir external lines
Line(13) = {10,11}; 
Line(14) = {11,12}; 
Line(15) = {12,13}; 
Line(16) = {13,10}; 
//inclided lines connecting reservoir to wellbore
Line(17) = {3,10};
Line(18) = {5,11};
Line(19) = {7,12};
Line(20) = {9,13};

//wellbore surfaces
//inner quad
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
// outer cylinder
Curve Loop(2) = {9, 10, 11, 12}; 

//4 loops
Curve Loop(3) = {1, 6, -9, -5};
Curve Loop(4) = {2, 7, -10, -6};
Curve Loop(5) = {3, 8, -11, -7};
Curve Loop(6) = {4, 5, -12, -8};
Plane Surface(2) = {3};
Plane Surface(3) = {4};
Plane Surface(4) = {5};
Plane Surface(5) = {6};

//reservoir surfaces
Curve Loop(7) = {17,13,-18,11};
Curve Loop(8) = {18,14,-19,12};
Curve Loop(9) = {19,15,-20,9};
Curve Loop(10)= {20,16,-17,10};

Plane Surface(6) = {7};
Plane Surface(7) = {8};
Plane Surface(8) = {9};
Plane Surface(9) = {10};

Transfinite Line {1, 2, 3, 4} = t_div;
Transfinite Line {9, 10, 11, 12} = t_div;
Transfinite Line {14, 16} = t_div; //reservoir farfield
Transfinite Line {13, 15} = t_div Using Bump b_res; //reservoir cap rock (top and bottom lines)
Transfinite Line {5, 6, 7, 8} = r_div; //radial lines of the wellbore
Transfinite Line {17, 18, 19, 20} = i_div Using Progression p_res; //lines connecting reservoir to wellbore

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

//extrude along x axis (wellbore cylinder)
extruded_1_ids() = Extrude {Lw, 0, 0} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Layers {h_div}; Recombine;
};
//Printf("%f",extruded_1_ids[12]);
//Transfinite Line {21,22,24,26,29,31,34,37} = h_div Using Bump 0.4; //tentative of making bump in wellbore (not working)
//Transfinite Line {41,43,46,49} = h_div;

//extrude along x (reservoir close to wellbore heel)
extruded_2_ids() = Extrude {-(Lr-Lw)/2, 0, 0} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Layers {h_div}; Recombine;
};

//extrude along x (reservoir close to wellbore toe)
extruded_3_ids() = Extrude {(Lr-Lw)/2, 0, 0} {
	Surface{extruded_1_ids[0]}; 
	Surface{extruded_1_ids[6]};
	Surface{extruded_1_ids[12]};
	Surface{extruded_1_ids[18]};
	Surface{extruded_1_ids[24]};
	Surface{extruded_1_ids[30]};
	Surface{extruded_1_ids[36]};
	Surface{extruded_1_ids[42]};
	Surface{extruded_1_ids[48]};
	Layers {h_div}; Recombine;
};

//remove the volume realted to the wellbore
Delete{ Point{1}; 
		Volume{extruded_1_ids[1]}; 
		Volume{extruded_1_ids[7]}; 
		Volume{extruded_1_ids[13]}; 
		Volume{extruded_1_ids[19]}; 
		Volume{extruded_1_ids[25]}; 
}

//set physical entites
Physical Curve("curve_wellbore",100) = {29};//just one line of the wellbore
Physical Curve("curve_toe",101) = {32,36,39,40};
Physical Curve("curve_heel",102) = {9,10,11,12};
Physical Surface("surface_wellbore_cylinder",103) = {16,20,23,25};//wellbore cylinder
Physical Surface("surface_wellbore_toe",104) = {14,18,21,24,26};//wellbore toe
Physical Surface("surface_wellbore_heel",105) = {1,2,3,4,5};//wellbore heel
Physical Surface("surface_farfield",106) = {66,37,95,31,60,89,59,62,65,67,94,96,88,91,43,47,50,53,55,72,76,82,84,79};//reservoir farfield
Physical Surface("cap_rock",107) = {57,28,86,63,34,92};//cap rock
all_volumes[] = Volume "*";
Physical Volume("volume_reservoir",108) = all_volumes[]; //reservoir

