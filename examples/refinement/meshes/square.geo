// Gmsh project created on Sun Feb 24 17:59:27 2019
SetFactory("OpenCASCADE");
Mesh.Algorithm = 8;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.5, -0.1, 0, 1.0};
//+
Point(6) = {1, 0.5, 0, 1.0};
//+
Point(7) = {0.5, 1.1, 0, 1.0};
//+
Point(8) = {0.1, 0.5, 0, 1.0};
//+
Spline(1) = {1, 5, 2};
//+
Spline(2) = {2, 6, 3};
//+
Spline(3) = {3, 7, 4};
//+
Spline(4) = {4, 8, 1};
//+
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Recombine Surface{5};
Physical Surface(1) = {5};
