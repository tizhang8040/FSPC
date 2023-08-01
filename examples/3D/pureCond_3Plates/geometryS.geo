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
NInt_FeBot = 50;
NExt_FeBot = 20;
Nz_FeBot   = 3;
clInt_FeBot = Lx_FeBot/NInt_FeBot;
clExt_FeBot = Lx_FeBot/NExt_FeBot;

NFSI_FeBot = Round(NInt_FeBot-Ly_FSI*(NInt_FeBot-NExt_FeBot)/Ly_FeBot);
//Printf("NFSI_FeBot = %f", NFSI_FeBot);
clFSI_FeBot = Lx_FeBot/NFSI_FeBot;



// Fe-bottom
// points
Point(1) = {X0-Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(2) = {X0+Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(3) = {X0+Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(4) = {X0-Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(5) = {X0+Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(6) = {X0-Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};

Point(7)  = {X0-Lx_FeBot/2+gxFSI,        Y0-Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(8)  = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0-Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(9)  = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0,          Z0_FeBot, clFSI_FeBot};
Point(10) = {X0-Lx_FeBot/2+gxFSI,        Y0,          Z0_FeBot, clFSI_FeBot};
Point(11) = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0+Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(12) = {X0-Lx_FeBot/2+gxFSI,        Y0+Ly_FSI/2, Z0_FeBot, clFSI_FeBot};

// lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {3, 5};
Line(6) = {5, 6};
Line(7) = {6, 4};

// curves
Curve Loop(1) = { 1, 2, 3, 4};
Curve Loop(2) = {-3, 5, 6, 7};

// surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// 2D mesh
//Transfinite Curve{1} = 40;
//Transfinite Curve{2} = 30;
//Transfinite Curve{3} = 40;
//Transfinite Curve{4} = 30;

//Transfinite Surface{1};
Recombine Surface{1, 2};

// volume
FeBot[] = Extrude {0., 0., Lz_FeBot} {
  Surface{1, 2}; Layers{Nz_FeBot}; Recombine;
};

// FSI
// points
Point(1) = {X0-Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(2) = {X0+Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(3) = {X0+Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(4) = {X0-Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(5) = {X0+Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(6) = {X0-Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};

Mesh 3;


// Groups ===============
//Physical Surface("Top") = {vol[0]};
//Physical Surface("FSInterface") = {1};
//Physical Volume("Solid") = {vol[1]};

//+
Recursive Delete {
  Point{26}; Curve{9}; Point{8}; Curve{32}; Curve{33}; Point{7}; Point{16}; Curve{34}; Point{12}; Curve{10}; Curve{12}; 
}
//+
Recursive Delete {
  Surface{1}; 
}
