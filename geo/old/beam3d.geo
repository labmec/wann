// Define the points of the polygon (cross section of the beam)
Point(1) = {0, 0, 0, 5.0};
Point(2) = {0, 30, 0, 5.0};
Point(3) = {0, 45, 20, 5.0};
Point(4) = {0, 60, 20, 5.0};
Point(5) = {0, 60, 40, 20.0};
Point(6) = {0, 30, 40, 5.0};
Point(7) = {0, 15, 20, 5.0};
Point(8) = {0, 0, 20, 20.0};

// Define the lines of the polygon
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// Create a line loop and a plane surface
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Extrude the surface by 150 units in the z-direction
Extrude {-150, 0, 0} {
    Surface{1};
}

// Define a physical group for the volume
Physical Volume("dom") = {1};

// Define a physical group for the boundary with plane surface 50
Physical Surface("plane") = {50};

// Define a physical group for point 5
Physical Point("axialforce") = {5};

// Define a physical group for point 8
Physical Point("2p") = {8};

// Set mesh size for all points
// Characteristic Length {:} = 5.0;
Characteristic Length {5,8} = 20.0;

// Force the mesh to be as uniform as possible
// Mesh.Algorithm = 1; // Use the Delaunay algorithm for 2D meshing
// Mesh.Algorithm3D = 1; // Use the Delaunay algorithm for 3D meshing
// Mesh.Optimize = 1; // Optimize the mesh
// Mesh.OptimizeNetgen = 1; // Use Netgen optimizer
// Mesh.ElementOrder = 1; // Use first-order elements
