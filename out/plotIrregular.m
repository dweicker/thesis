load u;
load x;
load y;

gridDelaunay = delaunay(x,y);
trimesh(gridDelaunay,x,y,u);