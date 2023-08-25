// parameters
Geometry.AutoCoherence = 0;	// keep duplicate elementary entities
// geometry in [mm]
gX0_RT = 20.;
// rigid tool
rExt_RT = 10.;			// radius of the rigid tool
// Fe-top
Lx_FeTop  = 110.;
Ly_FeTop  = 80.;
facLy1_FeTop = 1.0;		// factor for the Ly1 length
Ly1_FeTop = 2*rExt_RT*facLy1_FeTop;
Lz_FeTop  = 1.5;
// Al-middle (solid)
Lx_AlS = Lx_FeTop;
Ly_AlS = Ly_FeTop;
Lz_AlS = 3.;
// Al-middle (fluid)
Lx_AlF = 40.;
Ly_AlF = 30.;
// origin
X0 = 0.;
Y0 = 0.;
// Al-middle (solid)
Z0_AlS = 0.0;
// Fe-top
Z0_FeTop = Z0_AlS+Lz_AlS;
//caracteristic lengths
// Al-middle (solid)
NxInt_AlS = 50;
NxExt_AlS = 30;
Nz_AlS   = 3;
clInt_AlS = Lx_AlS/(NxInt_AlS-1);
clExt_AlS = Lx_AlS/(NxExt_AlS-1);
clFSI_AlS = clInt_AlS-Ly_AlF*(clInt_AlS-clExt_AlS)/Ly_AlS;
// Fe-top
//NxInt_FeTop = 80;
Ndy1_FeTop = 14;  //14
If (Ndy1_FeTop%2!=0)
    Ndy1_FeTop = Ndy1_FeTop + 1;
EndIf
clInt_FeTop = Ly1_FeTop/(Ndy1_FeTop-1);
Nx1_FeTop = Round(gX0_RT/clInt_FeTop);
Nx2_FeTop = Round(Lx_AlF/clInt_FeTop);
Nx3_FeTop = Round((Lx_FeTop-gX0_RT-Lx_AlF)/clInt_FeTop);
Printf("NxInt_FeTop = %f", (Nx1_FeTop+Nx2_FeTop+Nx3_FeTop));
NxExt_FeTop = 40;
Nz_FeTop = 3;
clExt_FeTop = Lx_FeTop/(NxExt_FeTop-1);
clFSI_FeTop = clInt_FeTop-Ly_AlF*(clInt_FeTop-clExt_FeTop)/Ly_FeTop;
Printf("clInt_FeTop = %f", clInt_FeTop);
Printf("clFSI_FeTop = %f", clFSI_FeTop);



// Al-middle (solid)
// points
Point(2001) = {X0-Lx_AlS/2, Y0-Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2002) = {X0+Lx_AlS/2, Y0-Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2003) = {X0-Lx_AlS/2, Y0,          Z0_AlS, clInt_AlS};
Point(2004) = {X0+Lx_AlS/2, Y0,          Z0_AlS, clInt_AlS};
Point(2005) = {X0-Lx_AlS/2, Y0+Ly_AlS/2, Z0_AlS, clExt_AlS};
Point(2006) = {X0+Lx_AlS/2, Y0+Ly_AlS/2, Z0_AlS, clExt_AlS};

Point(2007) = {X0-Lx_AlS/2+gX0_RT,        Y0-Ly_AlF/2, Z0_AlS, clFSI_AlS};
Point(2008) = {X0-Lx_AlS/2+gX0_RT+Lx_AlF, Y0-Ly_AlF/2, Z0_AlS, clFSI_AlS};
Point(2009) = {X0-Lx_AlS/2+gX0_RT,        Y0,          Z0_AlS, clInt_AlS};
Point(2010) = {X0-Lx_AlS/2+gX0_RT+Lx_AlF, Y0,          Z0_AlS, clInt_AlS};
Point(2011) = {X0-Lx_AlS/2+gX0_RT,        Y0+Ly_AlF/2, Z0_AlS, clFSI_AlS};
Point(2012) = {X0-Lx_AlS/2+gX0_RT+Lx_AlF, Y0+Ly_AlF/2, Z0_AlS, clFSI_AlS};

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

Point(4009) = {X0-Lx_FeTop/2+gX0_RT,        Y0-Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4010) = {X0-Lx_FeTop/2+gX0_RT+Lx_AlF, Y0-Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4011) = {X0-Lx_FeTop/2+gX0_RT,        Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4012) = {X0-Lx_FeTop/2+gX0_RT+Lx_AlF, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4013) = {X0-Lx_FeTop/2+gX0_RT,        Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4014) = {X0-Lx_FeTop/2+gX0_RT+Lx_AlF, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4015) = {X0-Lx_FeTop/2+gX0_RT,        Y0+Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4016) = {X0-Lx_FeTop/2+gX0_RT+Lx_AlF, Y0+Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};

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
Transfinite Curve{4017} = Ndy1_FeTop;
Transfinite Curve{4018} = Ndy1_FeTop;
Transfinite Curve{4019} = Ndy1_FeTop;
Transfinite Curve{4020} = Ndy1_FeTop;

Transfinite Curve{4007} = Nx1_FeTop+1;
Transfinite Curve{4015} = Nx1_FeTop+1;
Transfinite Curve{4021} = Nx2_FeTop+1;
Transfinite Curve{4022} = Nx2_FeTop+1;
Transfinite Curve{4003} = Nx3_FeTop+1;
Transfinite Curve{4011} = Nx3_FeTop+1;

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
Physical Surface("FSInterface") = {AlS_1[5],AlS_1[6],AlS_1[7],AlS_2[5],AlS_2[6],AlS_2[7], 4005,4006,4007};

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

Physical Volume("AlS") = {AlS_1[1],AlS_2[1]};
Physical Volume("FeTop") = {FeTop_1[1],FeTop_2[1],FeTop_3[1],FeTop_4[1],FeTop_5[1],FeTop_6[1],FeTop_7[1]};

