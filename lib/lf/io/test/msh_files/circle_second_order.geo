// Gmsh project created on Mon Apr 22 19:52:11 2019
//SetFactory("OpenCASCADE");

Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {0,1,0};
Point(4) = {-1,0,0};
Point(5) = {0,-1,0};

Circle(1) = {5,1,3};
Circle(2) = {3,1,5};


Line Loop(1) = {2, 1};
Plane Surface(1) = {1};

Physical Line(1) = {1, 2};
//+
Physical Surface(2) = {1};
