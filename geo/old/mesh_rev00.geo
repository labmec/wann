SetFactory("OpenCASCADE");

//geometry
Lw = 100; //[m] wellbore length 
Lr = 1000;//[m] reservoir width
nr = 10;//35;  //[ ] integer 
nw = 5;//50;  //[ ] integer
//mesh resolution
hw = Lw/20;
hr = Lr/10;

Point(1) = {0, 0, 0, hw}; //heel
Point(2) = {Lw, 0, 0, hw}; //toe

Point(3) = {-Lr,   -Lr, 0, hr};
Point(4) = {Lr+Lw, -Lr, 0, hr};
Point(5) = {Lr+Lw,  Lr, 0, hr};
Point(6) = {-Lr,    Lr, 0, hr};

Line(1) = {1, 2}; //wellbore
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 3}; 

Line(6) = {3, 1};
Line(7) = {4, 2};
Line(8) = {5, 2};
Line(9) = {6, 1}; 

Curve Loop(1) = {1,-7,-2,6};
Curve Loop(2) = {8,-7,3};
Curve Loop(3) = {8,-1,-9,-4};
Curve Loop(4) = {6,-9,5};

Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};
Surface(4) = {4};

Transfinite Line {1}  = nw Using Bump 0.1;
Transfinite Line {-6} = nr Using Progression 1.2;
Transfinite Line {-7} = nr Using Progression 1.2;
Transfinite Line {-8} = nr Using Progression 1.2;
Transfinite Line {-9} = nr Using Progression 1.2; 
//Transfinite Surface {1};

Physical Curve("far_field", 100) = {2,3,4,5};
Physical Curve("wellbore",101) = {1};
Physical Point("heel", 100) = {1};

Coherence Mesh;

algo_mesh_2D =  6;  //1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh.OptimizeThreshold = 0.5;
Mesh.Smoothing = 20;
Mesh 2;
OptimizeMesh "Gmsh";