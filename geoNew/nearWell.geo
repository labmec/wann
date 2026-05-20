Include "params.geo";

// Near-well region dimensions (in function of wellbore dimensions)
Lnw = Lw + Hr;
Wnw = Hr;
s = Sin(Pi/4.);

theta_div = h_div + 1; // Azimuthal division (of each circle quarter)

// Heel points
p1 = newp; Point(p1) = {0,     0,     0, 1.0}; // Center of the domain (wellbore heel)
p2 = newp; Point(p2) = {0, -Rw*s,  Rw*s, 1.0}; // Well
p3 = newp; Point(p3) = {0,  Rw*s,  Rw*s, 1.0}; // Well
p4 = newp; Point(p4) = {0,  Rw*s, -Rw*s, 1.0}; // Well
p5 = newp; Point(p5) = {0, -Rw*s, -Rw*s, 1.0}; // Well

// Circle arcs for heel
Circle(1) = {p2, p1, p3};
Circle(2) = {p3, p1, p4};
Circle(3) = {p4, p1, p5};
Circle(4) = {p5, p1, p2};

// Extrude to generate wellbore cylinder and toe
ids() = Extrude {Lw, 0, 0}{
   Line{1:4}; 
};

// --- Begin: Points and lines of the near-well region ---

// Close to wellbore heel
p6 = newp; Point(p6) = {-(Lnw-Lw)/2, -Wnw/2, Hr/2, 1.0};
p7 = newp; Point(p7) = {-(Lnw-Lw)/2, Wnw/2, Hr/2, 1.0}; 
p8 = newp; Point(p8) = {-(Lnw-Lw)/2, Wnw/2, -Hr/2, 1.0};
p9 = newp; Point(p9) = {-(Lnw-Lw)/2, -Wnw/2, -Hr/2, 1.0};

// Close to wellbore toe
p10 = newp; Point(p10) = {Lw+(Lnw-Lw)/2, -Wnw/2, Hr/2, 1.0};
p11 = newp; Point(p11) = {Lw+(Lnw-Lw)/2, Wnw/2, Hr/2, 1.0}; 
p12 = newp; Point(p12) = {Lw+(Lnw-Lw)/2, Wnw/2, -Hr/2, 1.0};
p13 = newp; Point(p13) = {Lw+(Lnw-Lw)/2, -Wnw/2, -Hr/2, 1.0};

// Lines aligned with the wellbore axis
l1 = newl; Line(l1) = {p6,p10}; 
l2 = newl; Line(l2) = {p7,p11}; 
l3 = newl; Line(l3) = {p8,p12}; 
l4 = newl; Line(l4) = {p9,p13}; 

// Lines close to the wellbore heel
l5 = newl; Line(l5) = {p6,p7}; 
l6 = newl; Line(l6) = {p7,p8}; 
l7 = newl; Line(l7) = {p8,p9}; 
l8 = newl; Line(l8) = {p9,p6}; 

// Lines close to the wellbore toe
l9 = newl; Line(l9)   = {p10,p11}; 
l10 = newl; Line(l10) = {p11,p12}; 
l11 = newl; Line(l11) = {p12,p13}; 
l12 = newl; Line(l12) = {p13,p10}; 

// Lines from the near-well region to the wellbore heel (diagonal lines)
l13 = newl; Line(l13) = {p2,p6};
l14 = newl; Line(l14) = {p3,p7};
l15 = newl; Line(l15) = {p4,p8};
l16 = newl; Line(l16) = {p5,p9};

// Lines from the near-well region to the wellbore toe (diagonal lines)
l17 = newl; Line(l17) = {p2+4,p10};
l18 = newl; Line(l18) = {p3+4,p11};
l19 = newl; Line(l19) = {p4+4,p12};
l20 = newl; Line(l20) = {p5+4,p13};

// --- End: Points and lines of the near-well region ---
// --- Begin: Curve loops and surfaces ---

// Surfaces of the near-well box
cl1 = newcl; Curve Loop(cl1) = {l5, l6, l7, l8};    // x-negative side, close to heel
cl2 = newcl; Curve Loop(cl2) = {l9, l10, l11, l12}; // x-positive side, close to toe
cl3 = newcl; Curve Loop(cl3) = {l1, l12, -l4, l8};  // y-negative side
cl4 = newcl; Curve Loop(cl4) = {l2, l10, -l3, l6};  // y-positive side
cl5 = newcl; Curve Loop(cl5) = {l5, l2, -l9, -l1};  // z-positive side, top
cl6 = newcl; Curve Loop(cl6) = {l7, l4, -l11, -l3}; // z-negative side, bottom
sf1 = news; Plane Surface(sf1) = {cl1};
sf2 = news; Plane Surface(sf2) = {cl2};
sf3 = news; Plane Surface(sf3) = {cl3};
sf4 = news; Plane Surface(sf4) = {cl4};
sf5 = news; Plane Surface(sf5) = {cl5};
sf6 = news; Plane Surface(sf6) = {cl6};

// Inclined surfaces close to the heel
cl7 = newcl; Curve Loop(cl7) = {1,-l14,-l5,l13};   // Top
cl8 = newcl; Curve Loop(cl8) = {2,-l15,-l6,l14};   // y-positive
cl9 = newcl; Curve Loop(cl9) = {3,-l16,-l7,l15};   // Bottom
cl10 = newcl; Curve Loop(cl10) = {4,-l13,-l8,l16}; // y-negative
sf7 = news; Surface(sf7) = {cl7};
sf8 = news; Surface(sf8) = {cl8};
sf9 = news; Surface(sf9) = {cl9};
sf10 = news; Surface(sf10) = {cl10};

// Inclined surfaces close to the toe
cl11 = newcl; Curve Loop(cl11) = {7,-l18,-l9,l17};   // Top
cl12 = newcl; Curve Loop(cl12) = {9,-l19,-l10,l18};  // y-positive
cl13 = newcl; Curve Loop(cl13) = {11,-l20,-l11,l19}; // Bottom
cl14 = newcl; Curve Loop(cl14) = {12,-l17,-l12,l20}; // y-negative
sf11 = news; Surface(sf11) = {cl11};
sf12 = news; Surface(sf12) = {cl12};
sf13 = news; Surface(sf13) = {cl13};
sf14 = news; Surface(sf14) = {cl14};

// Inclined surfaces (z-positive, top)
cl15 = newcl; Curve Loop(cl15) = {l13,5,-l17,-l1}; // y-negative
cl16 = newcl; Curve Loop(cl16) = {l14,6,-l18,-l2}; // y-positive
sf15 = news; Plane Surface(sf15) = {cl15};
sf16 = news; Plane Surface(sf16) = {cl16};

// Inclined surfaces (z-negative, bottom)
cl17 = newcl; Curve Loop(cl17) = {l16,10,-l20,-l4}; // y-negative
cl18 = newcl; Curve Loop(cl18) = {l15,8,-l19,-l3};  // y-positive
sf17 = news; Plane Surface(sf17) = {cl17};
sf18 = news; Plane Surface(sf18) = {cl18};

// Surfaces of the wellbore heel and wellbore toe
cl19 = newcl; Curve Loop(cl19) = {1,2,3,4}; // Heel
cl20 = newcl; Curve Loop(cl20) = {7,9,11,12}; // Toe
sf19 = news; Surface(sf19) = {cl19};
sf20 = news; Surface(sf20) = {cl20};

// --- End: Curve loops and surfaces ---
// --- Begin: Surface loops and Volumes ---

sfl1 = newsl; Surface Loop(sfl1) = {1,sf7,sf11,sf15,sf16,sf5};     // Top 
sfl2 = newsl; Surface Loop(sfl2) = {3,sf9,sf13,sf17,sf18,sf6};     // Bottom 
sfl3 = newsl; Surface Loop(sfl3) = {4,sf10,sf14,sf15,sf17,sf3};    // y-negative 
sfl4 = newsl; Surface Loop(sfl4) = {2,sf8,sf12,sf16,sf18,sf4};     // y-positive 
sfl5 = newsl; Surface Loop(sfl5) = {sf19,sf7,sf8,sf9,sf10,sf1};    // Wellbore heel 
sfl6 = newsl; Surface Loop(sfl6) = {sf20,sf11,sf12,sf13,sf14,sf2}; // Wellbore toe
vl1 = newv; Volume(vl1) = {sfl1};
vl2 = newv; Volume(vl2) = {sfl2};
vl3 = newv; Volume(vl3) = {sfl3};
vl4 = newv; Volume(vl4) = {sfl4};
vl5 = newv; Volume(vl5) = {sfl5};
vl6 = newv; Volume(vl6) = {sfl6};

// --- End: Surface loops and Volumes ---
// --- Begin: Meshing ---

// Apply transfinite progression to lines
Transfinite Line {l13,l14,l15,l16,l17,l18,l19,l20} = radial_div Using Progression p_res; // Diagonal lines of the near-well region
Transfinite Line {1, 2, 3, 4, 7, 9, 11, 12} = theta_div; // Circles of the heel and toe
Transfinite Line {l5,l6,l7,l8,l9,l10,l11,l12} = theta_div; // Lines of the near-well region (close to heel and toe)
Transfinite Line {l1,l2,l3,l4} = axial_div; // Horizontal near-well region lines (no bump)
Transfinite Line {5, 6, 8, 10} = axial_div Using Bump p_well; // Wellbore length (with bump)
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";

// Set physical entites for the final mesh
Physical Curve("curve_wellbore",300) = {5}; 
Physical Curve("curve_toe",301) = {7,9,11,12};
Physical Curve("curve_heel",302) = {1,2,3,4};
Physical Surface("surface_wellbore_cylinder",303) = {1,2,3,4}; 
Physical Surface("surface_wellbore_toe",304) = {sf20}; 
Physical Surface("surface_wellbore_heel",305) = {sf19};
Physical Point("point_heel",308) = {2};
Physical Point("point_toe",309) = {6};

// Set preliminary physical entities for the near-well region (to be merged with the reservoir mesh)
Physical Surface("surface_cap_rock_near_well",106) = {sf5,sf6};
Physical Volume("volume_near_well",110) = Volume "*";

Delete{ Point{p1}; } // Center point no longer needed

General.NumThreads = 4;
Mesh.Smoothing = 0;
Mesh 2;

// // Meshing smoothing loop
// For k In {0:40:1}
// 	OptimizeMesh "Laplace2D";
//    Printf("k = %f", k );
//    k++;
// EndFor

Mesh 3;
Coherence Mesh;
Save "nearWell.msh";
