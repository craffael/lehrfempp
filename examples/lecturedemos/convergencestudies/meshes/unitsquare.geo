Point(1) = {0.0, 0.0, 0, 1.0};
Point(2) = {1.0, 0.0, 0, 1.0};
Point(3) = {1.0, 1.0, 0, 1.0};
Point(4) = {0.0, 1.0, 0, 1.0};

Characteristic Length {4, 3, 2, 1} = 0.2;

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {3, 4, 1, 2};

Plane Surface(1) = {1};
