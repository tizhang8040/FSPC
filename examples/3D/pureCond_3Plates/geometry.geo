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
// Al-middle (solid)
Lx_AlS = Lx_FeTop+gExten;
Ly_AlS = Ly_FeTop+gExten;
Lz_AlS = 3.;
// Al-middle (fluid)
Lx_AlF = 40.;
Ly_AlF = 30.;
Lz_AlF = Lz_AlS;
// Fe-bottom
Lx_FeBot = Lx_AlS+gExten;
Ly_FeBot = Ly_AlS+gExten;
Lz_FeBot = 5.;
// origin
X0 = 0.;
Y0 = 0.;
// Fe-bottom
Z0_FeBot = 0.;
// Al-middle (solid)
Z0_AlS = Z0_FeBot+Lz_FeBot;
// Al-middle (fluid)
Z0_AlF = Z0_AlS;
// Fe-top
Z0_FeTop = Z0_AlS+Lz_AlS;
//caracteristic lengths
// Fe-bottom
NInt_FeBot = 30;
NExt_FeBot = 15;
Nz_FeBot   = 3;
clInt_FeBot = Lx_FeBot/NInt_FeBot;
clExt_FeBot = Lx_FeBot/NExt_FeBot;
clFSI_FeBot = clInt_FeBot-Ly_AlF*(clInt_FeBot-clExt_FeBot)/Ly_FeBot;
// Al-middle (solid)
NInt_AlS = 50;
NExt_AlS = 30;
Nz_AlS   = 3;
clInt_AlS = Lx_AlS/NInt_AlS;
clExt_AlS = Lx_AlS/NExt_AlS;
clFSI_AlS = clInt_AlS-Ly_AlF*(clInt_AlS-clExt_AlS)/Ly_AlS;
// Al-middle (fluid)
cl_AlF = 1.2;
// Fe-top
NInt_FeTop = 80;
NExt_FeTop = 50;
Nz_FeTop   = 3;
clInt_FeTop = Lx_FeTop/NInt_FeTop;
clExt_FeTop = Lx_FeTop/NExt_FeTop;

//NxFSI_FeTop = Lx_AlF/clInt_FeTop;
//If (NxFSI_FeTop%2==0)
//    NxFSI_FeTop = NxFSI_FeTop + 1;
//EndIf
//clInt_FeTop = Lx_AlF/NxFSI_FeTop;

clFSI_FeTop = clInt_FeTop-Ly_AlF*(clInt_FeTop-clExt_FeTop)/Ly_FeTop;
//Printf("clFSI_FeTop = %f", clFSI_FeTop);





// Fe-bottom
// points
Point(1001) = {X0-Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1002) = {X0+Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1003) = {X0-Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(1004) = {X0+Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(1005) = {X0-Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1006) = {X0+Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};

Point(1007) = {X0-Lx_FeBot/2+gxFSI,        Y0-Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1008) = {X0-Lx_FeBot/2+gxFSI+Lx_AlF, Y0-Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1009) = {X0-Lx_FeBot/2+gxFSI,        Y0,          Z0_FeBot, clInt_FeBot};
Point(1010) = {X0-Lx_FeBot/2+gxFSI+Lx_AlF, Y0,          Z0_FeBot, clInt_FeBot};
Point(1011) = {X0-Lx_FeBot/2+gxFSI,        Y0+Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1012) = {X0-Lx_FeBot/2+gxFSI+Lx_AlF, Y0+Ly_AlF/2, Z0_FeBot, clFSI_FeBot};

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


// Al-middle (solid)
// points
Point(2001) = {X0-Lx_AlS/2, Y0-Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2002) = {X0+Lx_AlS/2, Y0-Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2003) = {X0-Lx_AlS/2, Y0,          Z0_AlS, clInt_AlS};
Point(2004) = {X0+Lx_AlS/2, Y0,          Z0_AlS, clInt_AlS};
Point(2005) = {X0-Lx_AlS/2, Y0+Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2006) = {X0+Lx_AlS/2, Y0+Ly_AlS/2, Z0_AlS, clExt_AlS};

Point(2007) = {X0-Lx_AlS/2+gxFSI,        Y0-Ly_AlF/2, Z0_AlS, clFSI_AlS};
Point(2008) = {X0-Lx_AlS/2+gxFSI+Lx_AlF, Y0-Ly_AlF/2, Z0_AlS, clFSI_AlS};
Point(2009) = {X0-Lx_AlS/2+gxFSI,        Y0,          Z0_AlS, clInt_AlS};
Point(2010) = {X0-Lx_AlS/2+gxFSI+Lx_AlF, Y0,          Z0_AlS, clInt_AlS};
Point(2011) = {X0-Lx_AlS/2+gxFSI,        Y0+Ly_AlF/2, Z0_AlS, clFSI_AlS};
Point(2012) = {X0-Lx_AlS/2+gxFSI+Lx_AlF, Y0+Ly_AlF/2, Z0_AlS, clFSI_AlS};

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
AlS_1[] = Extrude {0., 0., Lz_AlS} {
  Surface{2001}; Layers{Nz_AlS}; Recombine;
};

AlS_2[] = Extrude {0., 0., Lz_AlS} {
  Surface{2002}; Layers{Nz_AlS}; Recombine;
};


// Al-middle (fluid)
// points
Point(3001) = {X0-Lx_AlS/2+gxFSI,        Y0-Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3002) = {X0-Lx_AlS/2+gxFSI+Lx_AlF, Y0-Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3003) = {X0-Lx_AlS/2+gxFSI+Lx_AlF, Y0+Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3004) = {X0-Lx_AlS/2+gxFSI,        Y0+Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3005) = {X0-Lx_AlS/2+gxFSI,        Y0-Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};
Point(3006) = {X0-Lx_AlS/2+gxFSI+Lx_AlF, Y0-Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};
Point(3007) = {X0-Lx_AlS/2+gxFSI+Lx_AlF, Y0+Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};
Point(3008) = {X0-Lx_AlS/2+gxFSI,        Y0+Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};

// lines
Line(3001) = {3001, 3002};
Line(3002) = {3002, 3003};
Line(3003) = {3003, 3004};
Line(3004) = {3004, 3001};

Line(3005) = {3005, 3006};
Line(3006) = {3006, 3007};
Line(3007) = {3007, 3008};
Line(3008) = {3008, 3005};

Line(3009) = {3001, 3005};
Line(3010) = {3002, 3006};
Line(3011) = {3003, 3007};
Line(3012) = {3004, 3008};

// curves
Curve Loop(3001) = {3001, 3002, 3003, 3004};
Curve Loop(3002) = {3005, 3006, 3007, 3008};
Curve Loop(3003) = {3001, 3010, -3005, -3009};
Curve Loop(3004) = {3002, 3011, -3006, -3010};
Curve Loop(3005) = {3003, 3012, -3007, -3011};
Curve Loop(3006) = {3004, 3009, -3008, -3012};

// surfaces
Plane Surface(3001) = {3001};
Plane Surface(3002) = {3002};
Plane Surface(3003) = {3003};
Plane Surface(3004) = {3004};
Plane Surface(3005) = {3005};
Plane Surface(3006) = {3006};

// volume
Surface Loop(3001) = {3001, 3002, 3003, 3004, 3005, 3006};
Volume(3001) = {3001};


// Fe-top
// points
Point(4001) = {X0-Lx_FeTop/2, Y0-Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(4002) = {X0+Lx_FeTop/2, Y0-Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(4003) = {X0-Lx_FeTop/2, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4004) = {X0+Lx_FeTop/2, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4005) = {X0-Lx_FeTop/2, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4006) = {X0+Lx_FeTop/2, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4007) = {X0-Lx_FeTop/2, Y0+Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(4008) = {X0+Lx_FeTop/2, Y0+Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};

Point(4009) = {X0-Lx_FeTop/2+gxFSI,        Y0-Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4010) = {X0-Lx_FeTop/2+gxFSI+Lx_AlF, Y0-Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4011) = {X0-Lx_FeTop/2+gxFSI,        Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4012) = {X0-Lx_FeTop/2+gxFSI+Lx_AlF, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4013) = {X0-Lx_FeTop/2+gxFSI,        Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4014) = {X0-Lx_FeTop/2+gxFSI+Lx_AlF, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4015) = {X0-Lx_FeTop/2+gxFSI,        Y0+Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4016) = {X0-Lx_FeTop/2+gxFSI+Lx_AlF, Y0+Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};

// lines
Line(4001) = {4001, 4002};
Line(4002) = {4002, 4004};
Line(4003) = {4004, 4012};
Line(4004) = {4012, 4010};
Line(4005) = {4010, 4009};
Line(4006) = {4009, 4011};
Line(4007) = {4011, 4003};
Line(4008) = {4003, 4001};

Line(4009) = {4007, 4008};
Line(4010) = {4008, 4006};
Line(4011) = {4006, 4014};
Line(4012) = {4014, 4016};
Line(4013) = {4016, 4015};
Line(4014) = {4015, 4013};
Line(4015) = {4013, 4005};
Line(4016) = {4005, 4007};

Line(4017) = {4003, 4005};
Line(4018) = {4011, 4013};
Line(4019) = {4012, 4014};
Line(4020) = {4004, 4006};
Line(4021) = {4011, 4012};
Line(4022) = {4013, 4014};

// curves
Curve Loop(4001) = {4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008};
Curve Loop(4002) = {4009, 4010, 4011, 4012, 4013, 4014, 4015, 4016};

Curve Loop(4003) = {4007, 4017, -4015, -4018};
Curve Loop(4004) = {4003, 4019, -4011, -4020};

Curve Loop(4005) = {4004, 4005, 4006, 4021};
Curve Loop(4006) = {4012, 4013, 4014, 4022};
Curve Loop(4007) = {4018, 4022, -4019, -4021};

// surfaces
Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};
Plane Surface(4004) = {4004};
Plane Surface(4005) = {4005};
Plane Surface(4006) = {4006};
Plane Surface(4007) = {4007};

// 2D mesh
//Transfinite Curve{1} = 40;
//Transfinite Curve{2} = 30;
//Transfinite Curve{3} = 40;
//Transfinite Curve{4} = 30;

Transfinite Surface{4003, 4004, 4007};
Recombine Surface{4001, 4002, 4003, 4004, 4005, 4006, 4007};

// volume
FeTop_1[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4001}; Layers{Nz_FeTop}; Recombine;
};

FeTop_2[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4002}; Layers{Nz_FeTop}; Recombine;
};

FeTop_3[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4003}; Layers{Nz_FeTop}; Recombine;
};

FeTop_4[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4004}; Layers{Nz_FeTop}; Recombine;
};

FeTop_5[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4005}; Layers{Nz_FeTop}; Recombine;
};

FeTop_6[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4006}; Layers{Nz_FeTop}; Recombine;
};

FeTop_7[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4007}; Layers{Nz_FeTop}; Recombine;
};


Mesh 3;


// Groups ===============
Physical Surface("InterfaceS") = {FeBot_3[0],FeBot_4[0], AlS_1[5],AlS_1[6],AlS_1[7],AlS_2[5],AlS_2[6],AlS_2[7], 4005,4006,4007};
Physical Surface("InterfaceF") = {3001, 3002, 3003, 3004, 3005, 3006};

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

Physical Surface("FeTop_Down") = {4001,4002,4003,4004};
Physical Surface("FeTop_Up") = {FeTop_1[0],FeTop_2[0],FeTop_3[0],FeTop_4[0],FeTop_5[0],FeTop_6[0],FeTop_7[0]};
Physical Surface("FeTop_Front") = {FeTop_1[9],FeTop_2[9],FeTop_3[3]};
Physical Surface("FeTop_End") = {FeTop_1[3],FeTop_2[3],FeTop_4[5]};
Physical Surface("FeTop_Sides") = {FeTop_1[2],FeTop_2[2]};

Physical Volume("FeBot") = {FeBot_1[1],FeBot_2[1],FeBot_3[1],FeBot_4[1]};
Physical Volume("AlS") = {AlS_1[1],AlS_2[1]};
Physical Volume("AlF") = {3001};
Physical Volume("FeTop") = {FeTop_1[1],FeTop_2[1],FeTop_3[1],FeTop_4[1],FeTop_5[1],FeTop_6[1],FeTop_7[1]};

