Mesh.RecombineAll = 1;
// grid size
lc = 1e-1;
// rectangle sizes
w = 1.0;
h = 1.0;

Point(1) = {0, 0, 0, lc};
Point(2) = {w, 0, 0, lc};
Point(3) = {w, h, 0, lc};
Point(4) = {0, h, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Mesh 2;
Save "rect_4.vtk";
