// parameters
L = 40 ;
H = 3 ;
Hs= 1.5 ; 
W = 30;

lc = 1.;

// points
Point(1) = {0, H, 0, lc};
Point(2) = {L, H, 0, lc};
Point(3) = {L, H, W, lc};
Point(4) = {0, H, W, lc};

// lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// curves
Curve Loop(1) = {1, 2, 3, 4};

// surfaces
Plane Surface(1) = {1};

// 2D mesh
//Transfinite Curve{1} = 40;
//Transfinite Curve{2} = 30;
//Transfinite Curve{3} = 40;
//Transfinite Curve{4} = 30;

//Transfinite Surface{1};
Recombine Surface{1};

// volume
vol[] = Extrude {0, Hs, 0} {
  Surface{1}; Layers{3}; Recombine;
};

Mesh 3;


// Groups ===============
Physical Surface("Top") = {vol[0]};
Physical Surface("FSInterface") = {1};
Physical Volume("Solid") = {vol[1]};

