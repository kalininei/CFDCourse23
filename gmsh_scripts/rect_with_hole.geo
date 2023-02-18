SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 1.0, 1.0, 0};

Rectangle(1);
// MeshSize{ PointsOf{ Rectangle{1}; } } = 0.01;

Mesh 2;
Save "rect.vtk";
