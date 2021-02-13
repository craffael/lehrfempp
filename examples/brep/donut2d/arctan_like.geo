// Gmsh project created on Sat Feb 13 14:36:47 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 10};
//+
Point(2) = {0.5, 0.85, 0, 10};
Point(3) = {0.8, 0.92, 0, 10};
//+
Point(4) = {1.7, 1, 0, 10};
Point(5) = {-1,-0.9,0,10};

Point(12) = {-0.5,-0.85,0,10};
Point(13) = {-0.8,-0.92,0,10};
Point(14) = {-1.7,-1,0,10};


//+
Spline(1) = {14,13,12,1, 2, 3,4};

Point(20) = {0.5,0.9,0,10};
Point(21) = {0,1,0,10};
//+
Line(2) = {4, 21};
//+
Line(3) = {21, 14};
//+
Curve Loop(1) = {3, 1, 2};
//+
Plane Surface(1) = {1};
