// Gmsh project created on Sun Sep 01 17:32:20 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1/Sqrt(2), 1/Sqrt(2), 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};

Circle(2) = {2, 1, 3};
//+
Line(3) = {3, 1};
//+
Curve Loop(1) = {1, 2, 3};
//+
Plane Surface(1) = {1};

Periodic Curve {1} = {-3} Rotate { { 0, 0, 1}, { 0, 0, 0}, -Pi/4};//+//+

Physical Point("origin") = {1};
Physical Curve("arc") = {2};
Physical Surface(3) = {1};
