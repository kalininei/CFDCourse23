// grid size
lc = 1e-1;
// rectangle sizes
w = 1.0;
h = 1.0;
z = 1.0;

Point(1) = {0, 0, 0, lc};
Point(2) = {w, 0, 0, lc};
Point(3) = {w, h, 0, lc};
Point(4) = {0, h, 0, lc};

Point(5) = {0, 0, z, lc};
Point(6) = {w, 0, z, lc};
Point(7) = {w, h, z, lc};
Point(8) = {0, h, z, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {1, 10, -5, -9};
Curve Loop(4) = {2, 11, -6, -10};
Curve Loop(5) = {3, 12, -7, -11};
Curve Loop(6) = {4, 9, -8, -12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Mesh 3;
Save "cube.vtk";
