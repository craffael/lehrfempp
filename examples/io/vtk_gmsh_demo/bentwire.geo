Point(1) = {0, 0, 0, 0.5};
Point(2) = {1, 0, 0, 0.5};
Point(3) = {0, 2, 0, 0.5};
Point(4) = {1, 2, 0, 0.5};
Point(5) = {2, 2, 0, 0.5};
Point(6) = {2, 3, 0, 0.5};
Point(7) = {2, 4, 0, 0.5};
Point(8) = {4, 4, 0, 0.5};
Point(9) = {4, 3, 0, 0.5};
Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {2, 4};
Line(4) = {6, 9};
Line(5) = {9, 8};
Line(6) = {8, 7};
Circle(7) = {4, 5, 6};
Circle(8) = {3, 5, 7};
Curve Loop(1) = {1, 2, 3, 7, 4, 5, 6, -8};
Plane Surface(1) = {1};
//+ Bottom contact
Physical Curve("Contact0") = {2};
//+ Right contact
Physical Curve("Contact1") = {5};
//+ Insulated boundary 
Physical Curve("Insulated") = {3, 7, 4, 6, 8, 1};
//+ Computational domain
Physical Surface("wire") = {1};
//+
Physical Point("bottomcontact") = {1, 2};
//+
Physical Point("leftcontact") = {8, 9};
