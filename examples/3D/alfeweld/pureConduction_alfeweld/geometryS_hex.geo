L  = 40;
WS = 5;
WC = 25;
H  = 3;
Hs = 1.5;

Ny = 2;

NxOut = 10;
NxInn = 20;
If (NxOut%2==0)
    NxOut = NxOut + 1;
EndIf
If (NxInn%2==0)
    NxInn = NxInn + 1;
EndIf
lcxOut = L/NxOut;
lcxInn = L/NxInn;

NzVar = 15;
NzInn = 10;
//Printf("%f", NzInn);
If (NzInn%2!=0)
    NzInn = NzInn + 1;
EndIf
//Printf("%f", NzInn);



// points ===============
Point(1) = {0, H, 0,        lcxOut};
Point(2) = {L, H, 0,        lcxOut};
Point(3) = {L, H, WS,       lcxInn};
Point(4) = {0, H, WS,       lcxInn};
Point(5) = {L, H, WS+WC,    lcxInn};
Point(6) = {0, H, WS+WC,    lcxInn};
Point(7) = {L, H, WS+WC+WS, lcxOut};
Point(8) = {0, H, WS+WC+WS, lcxOut};

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
//Mesh.Algorithm = 8; //11;  //"Frontal-Delaunay for quads" meshing algorithm (better 2D (planar) quadrilateral meshes)
//Mesh.RecombinationAlgorithm = 3;  //2 or 3; full-quad recombination algorithm to generate full-quad meshes
Mesh.RecombineAll = 1;
//RecombineMesh;
//Recombine Surface{1, 2, 3};
//Mesh.SubdivisionAlgorithm = 1;  //subdivide the resulting hybrid mesh to generate full-quad meshes

//Transfinite Curve{1} = NxOut Using Progression 1.;
//Transfinite Curve{2} = NzVar Using Progression proLcx;
//Transfinite Curve{3} = NxInn Using Progression 1.;
//Transfinite Curve{4} = NzVar Using Progression 1./proLcx;
//Transfinite Surface{1};
Transfinite Curve{3} = NxInn;
Transfinite Curve{5} = NzInn;
Transfinite Curve{6} = NxInn;
Transfinite Curve{7} = NzInn;
Transfinite Surface{2};

//Mesh 2;

// 3D mesh =================
FeTopS1[] = Extrude {0, Hs, 0} {
  Surface{1}; Layers{Ny}; Recombine;
};

FeTopC[] = Extrude {0, Hs, 0} {
  Surface{2}; Layers{Ny}; Recombine;
};

FeTopS2[] = Extrude {0, Hs, 0} {
  Surface{3}; Layers{Ny}; Recombine;
};

Mesh 3;

// Groups ===============
Physical Surface("FeTopUp") = {FeTopS1[0],FeTopC[0],FeTopS2[0]};
Physical Surface("FeTopLat") = {FeTopS1[2],FeTopS1[3],FeTopS1[5], FeTopC[3],FeTopC[5], FeTopS2[3],FeTopS2[4],FeTopS2[5]};
Physical Surface("FSInterface") = {1,2,3};
Physical Volume("FeTop") = {FeTopS1[1],FeTopC[1],FeTopS2[1]};

Save "geometryS_hex.msh";
