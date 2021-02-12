// Gmsh project created on Fri Feb 12 10:24:16 2021
//+

//+
SetFactory("OpenCASCADE");
Circle(1) = {-0, -0, 0, 1, 0, 2*Pi};
//+
Circle(2) = {-0, -0, 0, 0.9, 0, 2*Pi};

//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
