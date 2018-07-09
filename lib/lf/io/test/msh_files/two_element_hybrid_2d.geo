// Gmsh project created on Fri Jun 29 18:18:05 2018
SetFactory("OpenCASCADE");


Point(1) = {0, 0, 0, 10.0};
Point(2) = {1, 0, 0, 10.0};
Point(3) = {2,0,0,10.0};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};


Extrude{0,1,0} {
  Line{1};
  Layers{1};
  Recombine;
}

Line(6) = {3, 5};



Line Loop(2) = {4, -6, -2};
Plane Surface(2) = {2};


//+
Physical Point("physicalEntity1", 1) = {1};
//+
Physical Point("physicalEntity2", 2) = {1};

//+
Physical Surface("physicalEntity1", 1) = {2};
Physical Surface("physicalEntity3", 3) = {2};

Physical Line("diagonal") = {6};
Physical Surface("square") = {1};