SetFactory("OpenCASCADE");

//geometry
Lw = 400; //[m] wellbore length 
Lr = 1000;//[m] reservoir width

//mesh resolution
hw = Lw/40;
hr = Lr/16;

Point(1) = {0, 0, 0, hw}; //heel
Point(2) = {Lw, 0, 0, hw}; //toe

Line(1) = {1, 2}; //wellbore
Rectangle(1) = {-Lr, -Lr, 0, 2*Lr+Lw, 2*Lr, 0};

Line(6) = {3, 1};
Line(7) = {4, 2};
Line(8) = {5, 2};
Line(9) = {6, 1}; 

Line{1} In Surface{1};
Line{6} In Surface{1};
Line{7} In Surface{1};
Line{8} In Surface{1};
Line{9} In Surface{1};

MeshSize{3,4,5,6} = hr;

Physical Curve("far_field", 100) = {2,3,4,5};
Physical Curve("wellbore",101) = {1};
Physical Surface("dom",102) = {1};
Physical Point("heel", 100) = {1};
Physical Point("dedao", 105) = {2};

Coherence Mesh;

algo_mesh_2D =  1;  //1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay (default), 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
Mesh.OptimizeThreshold = 0.5;
Mesh.Smoothing = 3;
Mesh 2;
OptimizeMesh "Gmsh";
