// parameters
Geometry.AutoCoherence = 0;  // keep duplicate elementary entities
// geometry in [mm]
gExten = 0.;
gxFSI = 20.;
// Fe-top
Lx_FeTop  = 110.;
Ly_FeTop  = 80.;
Ly1_FeTop = 20.;
Lz_FeTop  = 1.5;
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
NInt_FeBot = 30;
NExt_FeBot = 15;
Nz_FeBot   = 3;
clInt_FeBot = Lx_FeBot/NInt_FeBot;
clExt_FeBot = Lx_FeBot/NExt_FeBot;
clFSI_FeBot = clInt_FeBot-Ly_FSI*(clInt_FeBot-clExt_FeBot)/Ly_FeBot;
//Printf("clFSI_FeBot = %f", clFSI_FeBot);
// Al-middle
NInt_Al = 50;
NExt_Al = 30;
Nz_Al   = 3;
clInt_Al = Lx_Al/NInt_Al;
clExt_Al = Lx_Al/NExt_Al;
clFSI_Al = clInt_Al-Ly_FSI*(clInt_Al-clExt_Al)/Ly_Al;
//Printf("clFSI_Al = %f", clFSI_Al);
// Fe-top
NInt_FeTop = 80;
NExt_FeTop = 50;
Nz_FeTop   = 3;
clInt_FeTop = Lx_FeTop/NInt_FeTop;
clExt_FeTop = Lx_FeTop/NExt_FeTop;

//NxFSI_FeTop = Lx_FSI/clInt_FeTop;
//If (NxFSI_FeTop%2==0)
//    NxFSI_FeTop = NxFSI_FeTop + 1;
//EndIf
//clInt_FeTop = Lx_FSI/NxFSI_FeTop;

clFSI_FeTop = clInt_FeTop-Ly_FSI*(clInt_FeTop-clExt_FeTop)/Ly_FeTop;
//Printf("clFSI_FeTop = %f", clFSI_FeTop);





// Fe-bottom
// points
Point(1001) = {X0-Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1002) = {X0+Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1003) = {X0-Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(1004) = {X0+Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(1005) = {X0-Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1006) = {X0+Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};

Point(1007) = {X0-Lx_FeBot/2+gxFSI,        Y0-Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(1008) = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0-Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(1009) = {X0-Lx_FeBot/2+gxFSI,        Y0,          Z0_FeBot, clInt_FeBot};
Point(1010) = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0,          Z0_FeBot, clInt_FeBot};
Point(1011) = {X0-Lx_FeBot/2+gxFSI,        Y0+Ly_FSI/2, Z0_FeBot, clFSI_FeBot};
Point(1012) = {X0-Lx_FeBot/2+gxFSI+Lx_FSI, Y0+Ly_FSI/2, Z0_FeBot, clFSI_FeBot};

// lines
Line(1001) = {1001, 1002};
Line(1002) = {1002, 1004};
Line(1003) = {1004, 1010};
Line(1004) = {1010, 1008};
Line(1005) = {1008, 1007};
Line(1006) = {1007, 1009};
Line(1007) = {1009, 1003};
Line(1008) = {1003, 1001};

Line(1009) = {1005, 1006};
Line(1010) = {1006, 1004};
Line(1011) = {1010, 1012};
Line(1012) = {1012, 1011};
Line(1013) = {1011, 1009};
Line(1014) = {1003, 1005};

Line(1015) = {1009, 1010};

// curves
Curve Loop(1001) = {1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008};
Curve Loop(1002) = {1009, 1010, 1003, 1011, 1012, 1013, 1007, 1014};
Curve Loop(1003) = {1004, 1005, 1006, 1015};
Curve Loop(1004) = {1011, 1012, 1013, 1015};

// surfaces
Plane Surface(1001) = {1001};
Plane Surface(1002) = {1002};
Plane Surface(1003) = {1003};
Plane Surface(1004) = {1004};

// 2D mesh
Recombine Surface{1001, 1002, 1003, 1004};

//Mesh 2;

// volume
FeBot_1[] = Extrude {0., 0., Lz_FeBot} {
  Surface{1001}; Layers{Nz_FeBot}; Recombine;
};

FeBot_2[] = Extrude {0., 0., Lz_FeBot} {
  Surface{1002}; Layers{Nz_FeBot}; Recombine;
};

FeBot_3[] = Extrude {0., 0., Lz_FeBot} {
  Surface{1003}; Layers{Nz_FeBot}; Recombine;
};

FeBot_4[] = Extrude {0., 0., Lz_FeBot} {
  Surface{1004}; Layers{Nz_FeBot}; Recombine;
};


// Al-middle
// points
Point(2001) = {X0-Lx_Al/2, Y0-Ly_Al/2, Z0_Al, clExt_Al};
Point(2002) = {X0+Lx_Al/2, Y0-Ly_Al/2, Z0_Al, clExt_Al};
Point(2003) = {X0-Lx_Al/2, Y0,         Z0_Al, clInt_Al};
Point(2004) = {X0+Lx_Al/2, Y0,         Z0_Al, clInt_Al};
Point(2005) = {X0-Lx_Al/2, Y0+Ly_Al/2, Z0_Al, clExt_Al};
Point(2006) = {X0+Lx_Al/2, Y0+Ly_Al/2, Z0_Al, clExt_Al};

Point(2007) = {X0-Lx_Al/2+gxFSI,        Y0-Ly_FSI/2, Z0_Al, clFSI_Al};
Point(2008) = {X0-Lx_Al/2+gxFSI+Lx_FSI, Y0-Ly_FSI/2, Z0_Al, clFSI_Al};
Point(2009) = {X0-Lx_Al/2+gxFSI,        Y0,          Z0_Al, clInt_Al};
Point(2010) = {X0-Lx_Al/2+gxFSI+Lx_FSI, Y0,          Z0_Al, clInt_Al};
Point(2011) = {X0-Lx_Al/2+gxFSI,        Y0+Ly_FSI/2, Z0_Al, clFSI_Al};
Point(2012) = {X0-Lx_Al/2+gxFSI+Lx_FSI, Y0+Ly_FSI/2, Z0_Al, clFSI_Al};

// lines
Line(2001) = {2001, 2002};
Line(2002) = {2002, 2004};
Line(2003) = {2004, 2010};
Line(2004) = {2010, 2008};
Line(2005) = {2008, 2007};
Line(2006) = {2007, 2009};
Line(2007) = {2009, 2003};
Line(2008) = {2003, 2001};

Line(2009) = {2005, 2006};
Line(2010) = {2006, 2004};
Line(2011) = {2010, 2012};
Line(2012) = {2012, 2011};
Line(2013) = {2011, 2009};
Line(2014) = {2003, 2005};

// curves
Curve Loop(2001) = {2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008};
Curve Loop(2002) = {2009, 2010, 2003, 2011, 2012, 2013, 2007, 2014};

// surfaces
Plane Surface(2001) = {2001};
Plane Surface(2002) = {2002};

// 2D mesh
Recombine Surface{2001, 2002};

// volume
AlS_1[] = Extrude {0., 0., Lz_Al} {
  Surface{2001}; Layers{Nz_Al}; Recombine;
};

AlS_2[] = Extrude {0., 0., Lz_Al} {
  Surface{2002}; Layers{Nz_Al}; Recombine;
};


// Fe-top
// points
Point(3001) = {X0-Lx_FeTop/2, Y0-Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(3002) = {X0+Lx_FeTop/2, Y0-Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(3003) = {X0-Lx_FeTop/2, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3004) = {X0+Lx_FeTop/2, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3005) = {X0-Lx_FeTop/2, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3006) = {X0+Lx_FeTop/2, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3007) = {X0-Lx_FeTop/2, Y0+Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(3008) = {X0+Lx_FeTop/2, Y0+Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};

Point(3009) = {X0-Lx_FeTop/2+gxFSI,        Y0-Ly_FSI/2,    Z0_FeTop, clFSI_FeTop};
Point(3010) = {X0-Lx_FeTop/2+gxFSI+Lx_FSI, Y0-Ly_FSI/2,    Z0_FeTop, clFSI_FeTop};
Point(3011) = {X0-Lx_FeTop/2+gxFSI,        Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3012) = {X0-Lx_FeTop/2+gxFSI+Lx_FSI, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3013) = {X0-Lx_FeTop/2+gxFSI,        Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3014) = {X0-Lx_FeTop/2+gxFSI+Lx_FSI, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(3015) = {X0-Lx_FeTop/2+gxFSI,        Y0+Ly_FSI/2,    Z0_FeTop, clFSI_FeTop};
Point(3016) = {X0-Lx_FeTop/2+gxFSI+Lx_FSI, Y0+Ly_FSI/2,    Z0_FeTop, clFSI_FeTop};

// lines
Line(3001) = {3001, 3002};
Line(3002) = {3002, 3004};
Line(3003) = {3004, 3012};
Line(3004) = {3012, 3010};
Line(3005) = {3010, 3009};
Line(3006) = {3009, 3011};
Line(3007) = {3011, 3003};
Line(3008) = {3003, 3001};

Line(3009) = {3007, 3008};
Line(3010) = {3008, 3006};
Line(3011) = {3006, 3014};
Line(3012) = {3014, 3016};
Line(3013) = {3016, 3015};
Line(3014) = {3015, 3013};
Line(3015) = {3013, 3005};
Line(3016) = {3005, 3007};

Line(3017) = {3003, 3005};
Line(3018) = {3011, 3013};
Line(3019) = {3012, 3014};
Line(3020) = {3004, 3006};
Line(3021) = {3011, 3012};
Line(3022) = {3013, 3014};

// curves
Curve Loop(3001) = {3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008};
Curve Loop(3002) = {3009, 3010, 3011, 3012, 3013, 3014, 3015, 3016};

Curve Loop(3003) = {3007, 3017, -3015, -3018};
Curve Loop(3004) = {3003, 3019, -3011, -3020};

Curve Loop(3005) = {3004, 3005, 3006, 3021};
Curve Loop(3006) = {3012, 3013, 3014, 3022};
Curve Loop(3007) = {3018, 3022, -3019, -3021};

// surfaces
Plane Surface(3001) = {3001};
Plane Surface(3002) = {3002};
Plane Surface(3003) = {3003};
Plane Surface(3004) = {3004};
Plane Surface(3005) = {3005};
Plane Surface(3006) = {3006};
Plane Surface(3007) = {3007};

// 2D mesh
//Transfinite Curve{1} = 40;
//Transfinite Curve{2} = 30;
//Transfinite Curve{3} = 40;
//Transfinite Curve{4} = 30;

Transfinite Surface{3003, 3004, 3007};
Recombine Surface{3001, 3002, 3003, 3004, 3005, 3006, 3007};

// volume
FeTop_1[] = Extrude {0., 0., Lz_FeTop} {
  Surface{3001}; Layers{Nz_FeTop}; Recombine;
};

FeTop_2[] = Extrude {0., 0., Lz_FeTop} {
  Surface{3002}; Layers{Nz_FeTop}; Recombine;
};

FeTop_3[] = Extrude {0., 0., Lz_FeTop} {
  Surface{3003}; Layers{Nz_FeTop}; Recombine;
};

FeTop_4[] = Extrude {0., 0., Lz_FeTop} {
  Surface{3004}; Layers{Nz_FeTop}; Recombine;
};

FeTop_5[] = Extrude {0., 0., Lz_FeTop} {
  Surface{3005}; Layers{Nz_FeTop}; Recombine;
};

FeTop_6[] = Extrude {0., 0., Lz_FeTop} {
  Surface{3006}; Layers{Nz_FeTop}; Recombine;
};

FeTop_7[] = Extrude {0., 0., Lz_FeTop} {
  Surface{3007}; Layers{Nz_FeTop}; Recombine;
};


Mesh 3;


// Groups ===============
Physical Surface("FSInterface") = {FeBot_3[0],FeBot_4[0], Al_1[5],Al_1[6],Al_1[7],Al_2[5],Al_2[6],Al_2[7], 3005,3006,3007};

Physical Surface("FeBot_Down") = {1001,1002,1003,1004};
Physical Surface("FeBot_Up") = {FeBot_1[0],FeBot_2[0]};
Physical Surface("FeBot_Front") = {FeBot_1[9],FeBot_2[9]};
Physical Surface("FeBot_End") = {FeBot_1[3],FeBot_2[3]};
Physical Surface("FeBot_Sides") = {FeBot_1[2],FeBot_2[2]};

Physical Surface("AlS_Down") = {2001,2002};
Physical Surface("AlS_Up") = {AlS_1[0],AlS_2[0]};
Physical Surface("AlS_Front") = {AlS_1[9],AlS_2[9]};
Physical Surface("AlS_End") = {AlS_1[3],AlS_2[3]};
Physical Surface("AlS_Sides") = {AlS_1[2],AlS_2[2]};

Physical Surface("FeTop_Down") = {3001,3002,3003,3004};
Physical Surface("FeTop_Up") = {FeTop_1[0],FeTop_2[0],FeTop_3[0],FeTop_4[0]};
Physical Surface("FeTop_Front") = {FeTop_1[9],FeTop_2[9],FeTop_3[3]};
Physical Surface("FeTop_End") = {FeTop_1[3],FeTop_2[3],FeTop_4[5]};
Physical Surface("FeTop_Sides") = {FeTop_1[2],FeTop_2[2]};

Physical Volume("FeBot") = {FeBot_1[1],FeBot_2[1],FeBot_3[1],FeBot_4[1]};
Physical Volume("Al") = {Al_1[1],Al_2[1]};
Physical Volume("FeTop") = {FeTop_1[1],FeTop_2[1],FeTop_3[1],FeTop_4[1],FeTop_5[1],FeTop_6[1],FeTop_7[1]};

