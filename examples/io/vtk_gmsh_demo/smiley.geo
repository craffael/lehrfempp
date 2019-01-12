// face
Point(1) = {0, 0, 0, 1.0};
Point(2) = {5, 0, 0, 1.0};
Circle(1) = {2, 1, 2};

// left eye
Point(3) = {-1.9, 1.9, 0, 1.0};
Point(4) = {-0.9, 1.9, 0, 1.0};
Circle(2) = {4,3,4};

// right eye
Point(5) = {1.9, 1.9, 0, 1.0};
Point(6) = {0.9, 1.9, 0, 1.0};
Circle(3) = {6,5,6};


// mouth
Point(7) = {-3, -1.5, 0, 1.0};
Point(8) = {3, -1.5, 0, 1.0};
Point(9) = {0, -3.6, 0, 1.0};
BSpline(4) = {7, 9, 8};

Line Loop(2) = {2};
Plane Surface(2) = {2};
Line Loop(3) = {3};
Plane Surface(3) = {3};
Line Loop(1) = {1};
Plane Surface(1) = {1,2,3};
Point {1} In Surface{1};
Line {4} In Surface {1};


Physical Point("nose", 1) = {1};
Physical Line("mouth", 5) = {4};
Physical Surface("face", 6) = {1,2,3};
Physical Surface("eyes", 7) = {2, 3};

