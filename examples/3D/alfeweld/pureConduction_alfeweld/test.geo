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
Point(1) = {0, H, 0,  100./30.};
Point(2) = {L, H, 0,  100./30.};
Point(3) = {L, H, WC, 100./120.};
Point(4) = {0, H, WC, 100./120.};

// lines ================
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// curves ===============
Curve Loop(1) = { 1, 2, 3, 4 };

// surfaces =============
Plane Surface(1) = {1};

// 2D mesh =================
Mesh.Algorithm = 8 //11;  //"Frontal-Delaunay for quads" meshing algorithm (better 2D (planar) quadrilateral meshes)
//Mesh.SubdivisionAlgorithm = 1;  //subdivide the resulting hybrid mesh to generate full-quad meshes
Mesh.RecombinationAlgorithm = 3;  //2 or 3; full-quad recombination algorithm to generate full-quad meshes
Mesh.RecombineAll = 1;


//Transfinite Curve{1} = 30;
Transfinite Curve{2} = 25 Using Progression 1/1.01;
//Transfinite Curve{3} = 120;
Transfinite Curve{4} = 25 Using Progression 1.01;
//Transfinite Surface{1};

//RecombineMesh;
//Recombine Surface{1, 2, 3};
//Recombine Surface{1,3};

Mesh 2;
