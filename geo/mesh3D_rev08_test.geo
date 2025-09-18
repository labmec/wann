SetFactory("OpenCASCADE");

Geometry.MatchMeshTolerance = 1e-12;
Geometry.Tolerance = 1e-16;
Mesh.ToleranceReferenceElement = 1e-10;
General.Terminal = 1;
Geometry.SnapPoints = 0;

//wellbore dimensions
Lw = 1; //[m] wellbore length //400 m
Rw = 0.2; //[m] wellbore radius
//Rw = 0.001/2; //[m] wellbore radius // FIXME this Rw value generates the error

//reservoir dimensions
Lr = 4; //[m] length //2000 m
Wr = 3; //[m] width  //1000 m
Hr = 3; //[m] height //1000 m

//wellbore eccentricity
ecc = 0; //[] 0<=ecc<=0.5 wellbore eccentricity (only over the z-axis) 
max_vert_disp = Hr/2-Rw; //[m]  

//wellbore division (chan
h_div = 8;  //horizontal division of the reservoir
r_div = 3;  //radial division (of each circle quarter)
l_div = 10; //axial division of the wellbore
p_res = 1.6;//progression coeficient of the reservoir mesh
p_well = 0.3;//progression coeficient of the wellbore mesh (bump scheme)

s = Sin(Pi/4.);

//points of the heel part {{{
p1 = newp; Point(p1) = {0,     0,     0 + ecc*max_vert_disp, 1.0};//center of the domain (wellbore heel)
p2 = newp; Point(p2) = {0, -Rw*s,  Rw*s + ecc*max_vert_disp, 1.0};//well
p3 = newp; Point(p3) = {0,  Rw*s,  Rw*s + ecc*max_vert_disp, 1.0};//well
p4 = newp; Point(p4) = {0,  Rw*s, -Rw*s + ecc*max_vert_disp, 1.0};//well
p5 = newp; Point(p5) = {0, -Rw*s, -Rw*s + ecc*max_vert_disp, 1.0};//well
//}}}
//circles and lines of the wellbore heel {{{
Circle(1) = {p2, p1, p3};
Circle(2) = {p3, p1, p4};
Circle(3) = {p4, p1, p5};
Circle(4) = {p5, p1, p2};
//}}}
//extrude to generate the wellbore toe (wellbore cylinder is also generated {{{
ids() = Extrude {Lw, 0, 0}{
   Line{1:4}; 
   //Layers {1};
   //Recombine;
};
//}}}
//create points and lines of the reservoir {{{
//close to wellbore heel
p6 = newp; Point(p6) = {-(Lr-Lw)/2, -Wr/2, Hr/2, 1.0};
p7 = newp; Point(p7) = {-(Lr-Lw)/2, Wr/2, Hr/2, 1.0}; 
p8 = newp; Point(p8) = {-(Lr-Lw)/2, Wr/2, -Hr/2, 1.0};
p9 = newp; Point(p9) = {-(Lr-Lw)/2, -Wr/2, -Hr/2, 1.0}; 
//close to wellbore toe
p10 = newp; Point(p10) = {Lw+(Lr-Lw)/2, -Wr/2, Hr/2, 1.0};
p11 = newp; Point(p11) = {Lw+(Lr-Lw)/2, Wr/2, Hr/2, 1.0}; 
p12 = newp; Point(p12) = {Lw+(Lr-Lw)/2, Wr/2, -Hr/2, 1.0};
p13 = newp; Point(p13) = {Lw+(Lr-Lw)/2, -Wr/2, -Hr/2, 1.0}; 
//lines aligned with the wellbore axis
l1 = newl; Line(l1) = {p6,p10}; 
l2 = newl; Line(l2) = {p7,p11}; 
l3 = newl; Line(l3) = {p8,p12}; 
l4 = newl; Line(l4) = {p9,p13}; 
//lines close to the wellbore heel
l5 = newl; Line(l5) = {p6,p7}; 
l6 = newl; Line(l6) = {p7,p8}; 
l7 = newl; Line(l7) = {p8,p9}; 
l8 = newl; Line(l8) = {p9,p6}; 
//lines close to the wellbore toe
l9 = newl; Line(l9)   = {p10,p11}; 
l10 = newl; Line(l10) = {p11,p12}; 
l11 = newl; Line(l11) = {p12,p13}; 
l12 = newl; Line(l12) = {p13,p10}; 
//lines from the reservoir to the wellbore heel (diagonal lines)
l13 = newl; Line(l13) = {p2,p6};
l14 = newl; Line(l14) = {p3,p7};
l15 = newl; Line(l15) = {p4,p8};
l16 = newl; Line(l16) = {p5,p9};
//lines from the reservoir to the wellbore toe (diagonal lines)
l17 = newl; Line(l17) = {p2+4,p10};
l18 = newl; Line(l18) = {p3+4,p11};
l19 = newl; Line(l19) = {p4+4,p12};
l20 = newl; Line(l20) = {p5+4,p13};
//}}}
//create curve loops and surfaces {{{
//res farfield
cl1 = newcl; Curve Loop(cl1) = {l5, l6, l7, l8}; //close to heel
cl2 = newcl; Curve Loop(cl2) = {l9, l10, l11, l12};//close to toe
cl3 = newcl; Curve Loop(cl3) = {l1, l12, -l4, l8};//y-negative side
cl4 = newcl; Curve Loop(cl4) = {l2, l10, -l3, l6};//y-positive side
sf1 = news; Plane Surface(sf1) = {cl1};
sf2 = news; Plane Surface(sf2) = {cl2};
sf3 = news; Plane Surface(sf3) = {cl3};
sf4 = news; Plane Surface(sf4) = {cl4};
//res cap rock
cl5 = newcl; Curve Loop(cl5) = {l5, l2, -l9, -l1};//top surface
cl6 = newcl; Curve Loop(cl6) = {l7, l4, -l11, -l3};//bottom surface
sf5 = news; Plane Surface(sf5) = {cl5};
sf6 = news; Plane Surface(sf6) = {cl6};
//inclined surfaces close to the heel
cl7 = newcl; Curve Loop(cl7) = {1,-l14,-l5,l13};//top
cl8 = newcl; Curve Loop(cl8) = {2,-l15,-l6,l14};//y-positive
cl9 = newcl; Curve Loop(cl9) = {3,-l16,-l7,l15};//bottom
cl10 = newcl; Curve Loop(cl10) = {4,-l13,-l8,l16};//y-negative
sf7 = news; Surface(sf7) = {cl7};
sf8 = news; Surface(sf8) = {cl8};
sf9 = news; Surface(sf9) = {cl9};
sf10 = news; Surface(sf10) = {cl10};
//inclined surfaces close to the toe
cl11 = newcl; Curve Loop(cl11) = {7,-l18,-l9,l17};//top
cl12 = newcl; Curve Loop(cl12) = {9,-l19,-l10,l18};//y-positive
cl13 = newcl; Curve Loop(cl13) = {11,-l20,-l11,l19};//bottom
cl14 = newcl; Curve Loop(cl14) = {12,-l17,-l12,l20};//y-negative
sf11 = news; Surface(sf11) = {cl11};
sf12 = news; Surface(sf12) = {cl12};
sf13 = news; Surface(sf13) = {cl13};
sf14 = news; Surface(sf14) = {cl14};
//inclined surfaces (z-positive, top)
cl15 = newcl;Curve Loop(cl15) = {l13,5,-l17,-l1};//y-negative
cl16 = newcl;Curve Loop(cl16) = {l14,6,-l18,-l2};//y-positive
sf15 = news; Plane Surface(sf15) = {cl15};
sf16 = news; Plane Surface(sf16) = {cl16};
//inclined surfaces (z-negative, bottom)
cl17 = newcl;Curve Loop(cl17) = {l16,10,-l20,-l4};//y-negative
cl18 = newcl;Curve Loop(cl18) = {l15,8,-l19,-l3};//y-positive
sf17 = news; Plane Surface(sf17) = {cl17};
sf18 = news; Plane Surface(sf18) = {cl18};
//surfaces of the wellbore heel and wellbore toe

//FIXME not sure why this curve loop command is changing the coordinates of point 2
cl19 = newcl; Curve Loop(cl19) = {1,2,3,4};//heel

//c[] = Point{p2};
//Printf("%f", c[1]); //to debug
//Printf("%f", c[2]); //to debug
//Abort;

cl20 = newcl; Curve Loop(cl20) = {7,9,11,12};//toe
sf19 = news; Surface(sf19) = {cl19};
sf20 = news; Surface(sf20) = {cl20};

//}}}


//create volumes {{{
sfl1 = newsl; Surface Loop(sfl1) = {1,sf7,sf11,sf15,sf16,sf5}; //top 
sfl2 = newsl; Surface Loop(sfl2) = {3,sf9,sf13,sf17,sf18,sf6}; //bottom 
sfl3 = newsl; Surface Loop(sfl3) = {4,sf10,sf14,sf15,sf17,sf3};//y negative 
sfl4 = newsl; Surface Loop(sfl4) = {2,sf8,sf12,sf16,sf18,sf4};//y positive 
sfl5 = newsl; Surface Loop(sfl5) = {sf19,sf7,sf8,sf9,sf10,sf1}; //wellbore heel 
sfl6 = newsl; Surface Loop(sfl6) = {sf20,sf11,sf12,sf13,sf14,sf2}; //wellbore heel 

vl1 = newv; Volume(vl1) = {sfl1};
vl2 = newv; Volume(vl2) = {sfl2};
vl3 = newv; Volume(vl3) = {sfl3};
vl4 = newv; Volume(vl4) = {sfl4};
vl5 = newv; Volume(vl5) = {sfl5};
vl6 = newv; Volume(vl6) = {sfl6};
///}}}


//Printf("%f", l1); //to debug

//apply transfinite progression to lines {{{
Transfinite Line {l13,l14,l15,l16,l17,l18,l19,l20} = h_div Using Progression p_res; //diagonal lines of the reservoir
Transfinite Line {1, 2, 3, 4, 7, 9, 11, 12} = r_div; //circles of the heel and toe 
Transfinite Line {l5,l6,l7,l8,l9,l10,l11,l12} = r_div; //lines of the reservoir (close to heel and toe)
Transfinite Line {l1,l2,l3,l4} = l_div; //horizotanl reservoir lines (no bump)
Transfinite Line {5, 6, 8, 10} = l_div Using Bump p_well; //wellbore length (with bump)
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
//}}}
////set physical entites {{{
Physical Curve("curve_wellbore",100) = {5};//just one line of the wellbore
Physical Curve("curve_toe",101) = {7,9,11,12};
Physical Curve("curve_heel",102) = {1,2,3,4};
Physical Surface("surface_wellbore_cylinder",103) = {1,2,3,4};//wellbore cylinder
Physical Surface("surface_wellbore_toe",104) = {sf20};//wellbore toe
Physical Surface("surface_wellbore_heel",105) = {sf19};//wellbore heel
Physical Surface("surface_farfield",106) = {sf1,sf2,sf3,sf4};//reservoir farfield
Physical Surface("surface_cap_rock",107) = {sf5,sf6};//cap rock
Physical Point("point_heel",108) = {2};
Physical Point("point_toe",109) = {6};
all_volumes[] = Volume "*";
Physical Volume("volume_reservoir",110) = all_volumes[]; //reservoir
//}}}

Delete{ Point{p1}; 
}

General.NumThreads = 4;
Mesh.Smoothing = 0;
Mesh 2;
For k In {0:40:1}
	OptimizeMesh "Laplace2D";
   Printf("k = %f", k );
   k++;
EndFor
Mesh 3;

Coherence Mesh;
Save "mesh3D_rev08_test.msh";
