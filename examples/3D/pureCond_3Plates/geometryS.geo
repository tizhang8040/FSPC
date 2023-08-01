// parameters
// geometry in [mm]
gExten = 5.;
gxFSI = 20.;
// Fe-top
Lx_FeTop = 110.;
Ly_FeTop = 60.;
Lz_FeTop = 1.5;
// Al-middle
Lx_Al = Lx_FeTop+gExten;
Ly_Al = Ly_FeTop+gExten;
Lz_Al = 3.;
// Fe-bottom
Lx_FeBot = Lx_Al+gExten;
Ly_FeBot = Ly_Al+gExten;
Lz_FeBot = 5.;
// FSI
Lx_FSI = 40.;
Ly_FSI = 30.;
Lz_FSI = Lz_Al;
// origin
X0 = 0.;
Y0 = 0.;
// Fe-bottom
Z0_FeBot = 0.;
// Al-middle
Z0_Al = Z0_FeBot+Lz_FeBot;
// Fe-top
Z0_FeTop = Z0_Al+Lz_Al;
//caracteristic lengths
// Fe-bottom
NInt_FeBot = 200;
NExt_FeBot = 10;
Nz_FeBot   = 3;
clInt_FeBot = Lx_FeBot/NInt_FeBot;
clExt_FeBot = Lx_FeBot/NExt_FeBot;

//NFSI_FeBot = Round(NInt_FeBot-Ly_FSI*(NInt_FeBot-NExt_FeBot)/Ly_FeBot);
//Printf("NFSI_FeBot = %f", NFSI_FeBot);
clFSI_FeBot = clInt_FeBot-Ly_FSI*(clInt_FeBot-clExt_FeBot)/Ly_FeBot;



// Fe-bottom
// points
Point(1) = {X0-Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(2) = {X0+Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(3) = {X0-Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(4) = {X0+Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(5) = {X0-Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(6) = {X0+Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};

Point(7)  = {X0-Lx_FeBot/2+gxFSI,        Y0-Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(8)  = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0-Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(9)  = {X0-Lx_FeBot/2+gxFSI,        Y0,          Z0_FeBot, clInt_FeBot};
Point(10) = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0,          Z0_FeBot, clInt_FeBot};
Point(11) = {X0-Lx_FeBot/2+gxFSI,        Y0+Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(12) = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0+Ly_FSI/2, Z0_FeBot, clFSI_FeBot};

// lines
Line(1) = {1,  2};
Line(2) = {2,  4};
Line(3) = {4,  10};
Line(4) = {10, 8};
Line(5) = {8,  7};
Line(6) = {7,  9};
Line(7) = {9,  3};
Line(8) = {3,  1};

Line(9)  = {5,  6};
Line(10) = {6,  4};
Line(11) = {10, 12};
Line(12) = {12, 11};
Line(13) = {11, 9};
Line(14) = {3,  5};

Line(15) = {9,  10};

// curves
Curve Loop(1) = { 1,  2,  3,  4,  5,   6,  7, 8};
Curve Loop(2) = { 9,  10, 3,  11, 12,  13, 7, 14};
Curve Loop(3) = { 4, 5, 6, 15};
Curve Loop(4) = { 11, 12, 13, 15};

// surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

// 2D mesh
//Transfinite Curve{1} = 40;
//Transfinite Curve{2} = 30;
//Transfinite Curve{3} = 40;
//Transfinite Curve{4} = 30;

//Transfinite Surface{1};
Recombine Surface{1, 2, 3, 4};

Mesh 2;

// volume
//FeBot[] = Extrude {0., 0., Lz_FeBot} {
//  Surface{1, 2}; Layers{Nz_FeBot}; Recombine;
//};


//Mesh 3;


// Groups ===============
//Physical Surface("Top") = {vol[0]};
//Physical Surface("FSInterface") = {1};
//Physical Volume("Solid") = {vol[1]};

