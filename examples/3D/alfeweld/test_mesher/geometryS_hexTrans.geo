L = 100;
WS = 20;
WC = 30;
H = 3;
Hs= 1.5; 

lcOut = 10.;
lcInn = 2.;

NzInn = WC/lcInn;
//Printf("%f", NzInn);
If (NzInn%2!=0)
    NzInn = NzInn + 1;
EndIf
//Printf("%f", NzInn);
lcInn = WC/NzInn;
NxInn = L/lcInn;
If (NxInn%2==0)
    NxInn = NxInn + 1;
EndIf

// points ===============
Point(1) = {0, H, 0,        lcOut};
Point(2) = {L, H, 0,        lcOut};
Point(3) = {L, H, WS,       lcInn};
Point(4) = {0, H, WS,       lcInn};
Point(5) = {L, H, WS+WC,    lcInn};
Point(6) = {0, H, WS+WC,    lcInn};
Point(7) = {L, H, WS+WC+WS, lcOut};
Point(8) = {0, H, WS+WC+WS, lcOut};


// lines ================
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};
Line(5)  = {3, 5};
Line(6)  = {5, 6};
Line(7)  = {6, 4};
Line(8)  = {5, 7};
Line(9)  = {7, 8};
Line(10) = {8, 6};

// curves ===============
Curve Loop(1) = { 1, 2, 3, 4 };
Curve Loop(2) = {-3, 5, 6, 7 };
Curve Loop(3) = {-6, 8, 9, 10};

// surfaces =============
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// 2D mesh =================
//Transfinite Curve{3} = NxInn;
//Transfinite Curve{5} = NzInn;
//Transfinite Curve{6} = NxInn;
//Transfinite Curve{7} = NzInn;

Transfinite Surface{2};


//RecombineMesh;
//Mesh.Algorithm = 8;  //"Frontal-Delaunay for quads" meshing algorithm (better 2D (planar) quadrilateral meshes)
//Mesh.SubdivisionAlgorithm = 1;  //subdivide the resulting hybrid mesh to generate full-quad meshes
//Mesh.RecombinationAlgorithm = 3;  //2 or 3; full-quad recombination algorithm to generate full-quad meshes


//Recombine Surface{1, 2, 3};
Recombine Surface{:};

Mesh 2;


