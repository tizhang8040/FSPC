// parameters
Geometry.AutoCoherence = 0;	// keep duplicate elementary entities
// geometry in [mm]
// rigid tool
rExt_RT = 12.0;			// radius of the rigid tool
// Al-middle (fluid)
Lx_AlF = 110.;
Ly_AlF = 20.;
gx_AlF = 12.5;			// gap in x direction between left edges of AlS-AlF
// Fe-top
Lx_FeTop  = 135.;
Ly_FeTop  = 60.;
Lz_FeTop  = 1.5;
// Al-middle (solid)
Lx_AlS = Lx_FeTop;
Ly_AlS = Ly_FeTop;
Lz_AlS = 3.;
// Fe-bottom
Lx_FeBot = Lx_AlS;
Ly_FeBot = Ly_AlS;
Lz_FeBot = 5.;
Lx_FeBot_Central = 122.5;
// origin
X0 = 0.;
Y0 = 0.;
// Fe-bottom
Z0_FeBot = 0.;
// Al-middle (solid)
Z0_AlS = Z0_FeBot+Lz_FeBot;
// Fe-top
Z0_FeTop = Z0_AlS+Lz_AlS;
//caracteristic lengths
// Fe-bottom
NxInt_FeBot = 40;
NxExt_FeBot = 20;
Nz_FeBot    = 2;
clInt_FeBot = Lx_FeBot/NxInt_FeBot;
clExt_FeBot = Lx_FeBot/NxExt_FeBot;
clFSI_FeBot = clInt_FeBot-Ly_AlF*(clInt_FeBot-clExt_FeBot)/Ly_FeBot;
// Al-middle (solid)
NxInt_AlS = 70;
NxExt_AlS = 30;
Nz_AlS   = 2;
clInt_AlS = Lx_AlS/NxInt_AlS;
clExt_AlS = Lx_AlS/NxExt_AlS;
clFSI_AlS = clInt_AlS-Ly_AlF*(clInt_AlS-clExt_AlS)/Ly_AlS;
clMid_AlS = clInt_AlS-gx_AlF*(clInt_AlS-clExt_AlS)/Ly_AlS;
clMid_AlS = Min(clMid_AlS,clExt_AlS);
// Fe-top
//NxInt_FeTop = 80;
Ny1_FeTop = 12;
If (Ny1_FeTop%2==0)
    Ny1_FeTop = Ny1_FeTop + 1;
EndIf
clInt_FeTop = Ly_AlF/Ny1_FeTop;
Nx2_FeTop = Round(Lx_AlF/clInt_FeTop);
//Printf("NxInt_FeTop = %f", (Nx1_FeTop+Nx2_FeTop+Nx3_FeTop));
NxExt_FeTop = 40;
Nz_FeTop = 2;
clExt_FeTop = Lx_FeTop/NxExt_FeTop;
clFSI_FeTop = clInt_FeTop-Ly_AlF*(clInt_FeTop-clExt_FeTop)/Ly_FeTop;
clMid_FeTop = clInt_FeTop-gx_AlF*(clInt_FeTop-clExt_FeTop)/Ly_FeTop;
clMid_FeTop = Min(clMid_FeTop,clExt_FeTop);
//Printf("clInt_FeTop = %f", clInt_FeTop);
//Printf("clFSI_FeTop = %f", clFSI_FeTop);



// Fe-bottom
// points
Point(1001) = {X0-Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1002) = {X0+Lx_FeBot/2, Y0-Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1003) = {X0-Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(1004) = {X0+Lx_FeBot/2, Y0,            Z0_FeBot, clInt_FeBot};
Point(1005) = {X0-Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};
Point(1006) = {X0+Lx_FeBot/2, Y0+Ly_FeBot/2, Z0_FeBot, clExt_FeBot};

Point(1007) = {X0-Lx_FeBot_Central/2,       Y0-Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1008) = {X0-Lx_FeBot/2+gx_AlF,        Y0-Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1009) = {X0-Lx_FeBot/2+gx_AlF+Lx_AlF, Y0-Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1010) = {X0+Lx_FeBot_Central/2,       Y0-Ly_AlF/2, Z0_FeBot, clFSI_FeBot};

Point(1011) = {X0-Lx_FeBot_Central/2,       Y0, Z0_FeBot, clInt_FeBot};
Point(1012) = {X0-Lx_FeBot/2+gx_AlF,        Y0, Z0_FeBot, clInt_FeBot};
Point(1013) = {X0-Lx_FeBot/2+gx_AlF+Lx_AlF, Y0, Z0_FeBot, clInt_FeBot};
Point(1014) = {X0+Lx_FeBot_Central/2,       Y0, Z0_FeBot, clInt_FeBot};

Point(1015) = {X0-Lx_FeBot_Central/2,       Y0+Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1016) = {X0-Lx_FeBot/2+gx_AlF,        Y0+Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1017) = {X0-Lx_FeBot/2+gx_AlF+Lx_AlF, Y0+Ly_AlF/2, Z0_FeBot, clFSI_FeBot};
Point(1018) = {X0+Lx_FeBot_Central/2,       Y0+Ly_AlF/2, Z0_FeBot, clFSI_FeBot};

// lines
Line(1001) = {1001, 1002};
Line(1002) = {1002, 1004};
Line(1003) = {1004, 1014};
Line(1004) = {1014, 1010};
Line(1005) = {1010, 1009};
Line(1006) = {1009, 1008};
Line(1007) = {1008, 1007};
Line(1008) = {1007, 1011};
Line(1009) = {1011, 1003};
Line(1010) = {1003, 1001};

Line(1011) = {1005, 1006};
Line(1012) = {1006, 1004};
Line(1014) = {1014, 1018};
Line(1015) = {1018, 1017};
Line(1016) = {1017, 1016};
Line(1017) = {1016, 1015};
Line(1018) = {1015, 1011};
Line(1020) = {1003, 1005};

Line(1021) = {1011, 1012};
Line(1022) = {1012, 1008};
Line(1023) = {1012, 1016};
Line(1024) = {1012, 1013};
Line(1025) = {1009, 1013};
Line(1026) = {1017, 1013};
Line(1027) = {1013, 1014};

// curves
Curve Loop(1001) = {1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010};
Curve Loop(1002) = {1011, 1012, 1003, 1014, 1015, 1016, 1017, 1018, 1009, 1020};
Curve Loop(1003) = {1007, 1008, 1021, 1022};
Curve Loop(1004) = {1017, 1018, 1021, 1023};
Curve Loop(1005) = {1006,-1022, 1024,-1025};
Curve Loop(1006) = {1016,-1023, 1024,-1026};
Curve Loop(1007) = {1005, 1025, 1027, 1004};
Curve Loop(1008) = {1015, 1026, 1027, 1014};

// surfaces
For k In {1001:1008:1}
   Plane Surface(k) = {k};
EndFor

// 2D mesh
Recombine Surface{1001,1002,1003,1004,1005,1006,1007,1008};
//Mesh 2;

// volume
FeBot_1[] = Extrude {0.,0.,Lz_FeBot} {Surface{1001}; Layers{Nz_FeBot}; Recombine;};
FeBot_2[] = Extrude {0.,0.,Lz_FeBot} {Surface{1002}; Layers{Nz_FeBot}; Recombine;};
FeBot_3[] = Extrude {0.,0.,Lz_FeBot} {Surface{1003}; Layers{Nz_FeBot}; Recombine;};
FeBot_4[] = Extrude {0.,0.,Lz_FeBot} {Surface{1004}; Layers{Nz_FeBot}; Recombine;};
FeBot_5[] = Extrude {0.,0.,Lz_FeBot} {Surface{1005}; Layers{Nz_FeBot}; Recombine;};
FeBot_6[] = Extrude {0.,0.,Lz_FeBot} {Surface{1006}; Layers{Nz_FeBot}; Recombine;};
FeBot_7[] = Extrude {0.,0.,Lz_FeBot} {Surface{1007}; Layers{Nz_FeBot}; Recombine;};
FeBot_8[] = Extrude {0.,0.,Lz_FeBot} {Surface{1008}; Layers{Nz_FeBot}; Recombine;};



// Al-middle (solid)
// points
Point(2001) = {X0-Lx_AlS/2, Y0-Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2002) = {X0+Lx_AlS/2, Y0-Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2003) = {X0-Lx_AlS/2, Y0,          Z0_AlS, clMid_AlS};
Point(2004) = {X0+Lx_AlS/2, Y0,          Z0_AlS, clMid_AlS};
Point(2005) = {X0-Lx_AlS/2, Y0+Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2006) = {X0+Lx_AlS/2, Y0+Ly_AlS/2, Z0_AlS, clExt_AlS};

Point(2007) = {X0-Lx_AlS/2+gx_AlF,        Y0-Ly_AlF/2, Z0_AlS, clInt_AlS};
Point(2008) = {X0-Lx_AlS/2+gx_AlF+Lx_AlF, Y0-Ly_AlF/2, Z0_AlS, clInt_AlS};
Point(2009) = {X0-Lx_AlS/2+gx_AlF,        Y0,          Z0_AlS, clInt_AlS};
Point(2010) = {X0-Lx_AlS/2+gx_AlF+Lx_AlF, Y0,          Z0_AlS, clInt_AlS};
Point(2011) = {X0-Lx_AlS/2+gx_AlF,        Y0+Ly_AlF/2, Z0_AlS, clInt_AlS};
Point(2012) = {X0-Lx_AlS/2+gx_AlF+Lx_AlF, Y0+Ly_AlF/2, Z0_AlS, clInt_AlS};

// lines
Line(2001) = {2001, 2002};
Line(2002) = {2002, 2004};
Line(2003) = {2004, 2010};
Line(2004) = {2010, 2008};
Line(2005) = {2008, 2007};
Line(2006) = {2007, 2009};
Line(2007) = {2009, 2003};
Line(2008) = {2003, 2001};

Line(2011) = {2005, 2006};
Line(2012) = {2006, 2004};
Line(2014) = {2010, 2012};
Line(2015) = {2012, 2011};
Line(2016) = {2011, 2009};
Line(2018) = {2003, 2005};

// curves
Curve Loop(2001) = {2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008};
Curve Loop(2002) = {2011, 2012, 2003, 2014, 2015, 2016, 2007, 2018};

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



// Fe-top
// points
Point(4001) = {X0-Lx_FeTop/2, Y0-Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(4002) = {X0+Lx_FeTop/2, Y0-Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(4003) = {X0-Lx_FeTop/2, Y0,             Z0_FeTop, clMid_FeTop};
Point(4004) = {X0+Lx_FeTop/2, Y0,             Z0_FeTop, clMid_FeTop};
Point(4005) = {X0-Lx_FeTop/2, Y0+Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};
Point(4006) = {X0+Lx_FeTop/2, Y0+Ly_FeTop/2,  Z0_FeTop, clExt_FeTop};

Point(4007) = {X0-Lx_FeTop/2+gx_AlF,        Y0-Ly_AlF/2,    Z0_FeTop, clInt_FeTop};
Point(4008) = {X0-Lx_FeTop/2+gx_AlF+Lx_AlF, Y0-Ly_AlF/2,    Z0_FeTop, clInt_FeTop};
Point(4009) = {X0-Lx_FeTop/2+gx_AlF,        Y0,             Z0_FeTop, clInt_FeTop};
Point(4010) = {X0-Lx_FeTop/2+gx_AlF+Lx_AlF, Y0,             Z0_FeTop, clInt_FeTop};
Point(4011) = {X0-Lx_FeTop/2+gx_AlF,        Y0+Ly_AlF/2,    Z0_FeTop, clInt_FeTop};
Point(4012) = {X0-Lx_FeTop/2+gx_AlF+Lx_AlF, Y0+Ly_AlF/2,    Z0_FeTop, clInt_FeTop};

// lines
Line(4001) = {4001, 4002};
Line(4002) = {4002, 4004};
Line(4003) = {4004, 4010};
Line(4004) = {4010, 4008};
Line(4005) = {4008, 4007};
Line(4006) = {4007, 4009};
Line(4007) = {4009, 4003};
Line(4008) = {4003, 4001};

Line(4011) = {4005, 4006};
Line(4012) = {4006, 4004};
Line(4014) = {4010, 4012};
Line(4015) = {4012, 4011};
Line(4016) = {4011, 4009};
Line(4018) = {4003, 4005};

Line(4021) = {4007, 4011};
Line(4022) = {4008, 4012};

// curves
Curve Loop(4001) = {4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008};
Curve Loop(4002) = {4011, 4012, 4003, 4014, 4015, 4016, 4007, 4018};

Curve Loop(4003) = {4005, 4021, -4015, -4022};

// surfaces
Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};

// 2D mesh
Transfinite Curve{4005} = Nx2_FeTop+1;
Transfinite Curve{4015} = Nx2_FeTop+1;

Transfinite Curve{4021} = Ny1_FeTop+1;
Transfinite Curve{4022} = Ny1_FeTop+1;

Transfinite Surface{4003};
Recombine Surface{4001, 4002, 4003};

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


Mesh 3;




// Groups
Physical Surface("FSInterface") = {FeBot_5[0],FeBot_6[0], AlS_1[5],AlS_1[6],AlS_1[7],AlS_2[5],AlS_2[6],AlS_2[7], 4003};

Physical Surface("FeBot_Down") = {1001,1002,1003,1004,1005,1006,1007,1008};
Physical Surface("FeBot_DownCentral") = {1003,1004,1005,1006,1007,1008};
Physical Surface("FeBot_Up") = {FeBot_1[0],FeBot_2[0],FeBot_3[0],FeBot_4[0],FeBot_5[0],FeBot_6[0],FeBot_7[0],FeBot_8[0]};
Physical Surface("FeBot_Front") = {FeBot_1[11],FeBot_2[11]};
Physical Surface("FeBot_End") = {FeBot_1[3],FeBot_2[3]};
Physical Surface("FeBot_Sides") = {FeBot_1[2],FeBot_2[2]};

Physical Surface("AlS_Down") = {2001,2002};
Physical Surface("AlS_Up") = {AlS_1[0],AlS_2[0]};
Physical Surface("AlS_Front") = {AlS_1[9],AlS_2[9]};
Physical Surface("AlS_End") = {AlS_1[3],AlS_2[3]};
Physical Surface("AlS_Sides") = {AlS_1[2],AlS_2[2]};

Physical Surface("FeTop_Down") = {4001,4002,4003};
Physical Surface("FeTop_Up") = {FeTop_1[0],FeTop_2[0],FeTop_3[0]};
Physical Surface("FeTop_UpCentral") = {FeTop_3[0]};
Physical Surface("FeTop_Front") = {FeTop_1[9],FeTop_2[9]};
Physical Surface("FeTop_End") = {FeTop_1[3],FeTop_2[3]};
Physical Surface("FeTop_Sides") = {FeTop_1[2],FeTop_2[2]};

Physical Volume("FeBot") = {FeBot_1[1],FeBot_2[1],FeBot_3[1],FeBot_4[1],FeBot_5[1],FeBot_6[1],FeBot_7[1],FeBot_8[1]};
Physical Volume("AlS") = {AlS_1[1],AlS_2[1]};
Physical Volume("FeTop") = {FeTop_1[1],FeTop_2[1],FeTop_3[1]};
