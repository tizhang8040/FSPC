// parameters
// geometry in [mm]
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
gx_AlF = 20.;			// gap in x direction between left edges of AlS-AlF
// origin
X0 = 0.;
Y0 = 0.;
// Fe-bottom
// Al-middle (solid)
Z0_AlS = 0.;
// Fe-top
Z0_FeTop = Z0_AlS+Lz_AlS;
//caracteristic lengths
// Fe-top
//NxInt_FeTop = 80;
Ny1_FeTop = 9;
If (Ny1_FeTop%2==0)
    Ny1_FeTop = Ny1_FeTop + 1;
EndIf
clInt_FeTop = Ly1_FeTop/Ny1_FeTop;
Nx1_FeTop = Round(gx_AlF/clInt_FeTop);
Nx2_FeTop = Round(Lx_AlF/clInt_FeTop);
Nx3_FeTop = Round((Lx_FeTop-gx_AlF-Lx_AlF)/clInt_FeTop);
Printf("NxInt_FeTop = %f", (Nx1_FeTop+Nx2_FeTop+Nx3_FeTop));
NxExt_FeTop = 25;
Nz_FeTop = 2;
clExt_FeTop = Lx_FeTop/NxExt_FeTop;
clFSI_FeTop = clInt_FeTop-Ly_AlF*(clInt_FeTop-clExt_FeTop)/Ly_FeTop;
Printf("clInt_FeTop = %f", clInt_FeTop);
Printf("clFSI_FeTop = %f", clFSI_FeTop);



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

Point(4009) = {X0-Lx_FeTop/2+gx_AlF,        Y0-Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4010) = {X0-Lx_FeTop/2+gx_AlF+Lx_AlF, Y0-Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4011) = {X0-Lx_FeTop/2+gx_AlF,        Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4012) = {X0-Lx_FeTop/2+gx_AlF+Lx_AlF, Y0-Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4013) = {X0-Lx_FeTop/2+gx_AlF,        Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4014) = {X0-Lx_FeTop/2+gx_AlF+Lx_AlF, Y0+Ly1_FeTop/2, Z0_FeTop, clInt_FeTop};
Point(4015) = {X0-Lx_FeTop/2+gx_AlF,        Y0+Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};
Point(4016) = {X0-Lx_FeTop/2+gx_AlF+Lx_AlF, Y0+Ly_AlF/2,    Z0_FeTop, clFSI_FeTop};

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
Transfinite Curve{4017} = Ny1_FeTop+1;
Transfinite Curve{4018} = Ny1_FeTop+1;
Transfinite Curve{4019} = Ny1_FeTop+1;
Transfinite Curve{4020} = Ny1_FeTop+1;

Transfinite Curve{4007} = Nx1_FeTop+1;
Transfinite Curve{4015} = Nx1_FeTop+1;
Transfinite Curve{4021} = Nx2_FeTop+1;
Transfinite Curve{4022} = Nx2_FeTop+1;
Transfinite Curve{4003} = Nx3_FeTop+1;
Transfinite Curve{4011} = Nx3_FeTop+1;

Transfinite Surface{4003, 4004, 4007};
//Recombine Surface{4001, 4002, 4003, 4004, 4005, 4006, 4007};

// volume
FeTop_1[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4001}; Layers{Nz_FeTop};
};

FeTop_2[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4002}; Layers{Nz_FeTop};
};

FeTop_3[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4003}; Layers{Nz_FeTop};
};

FeTop_4[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4004}; Layers{Nz_FeTop};
};

FeTop_5[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4005}; Layers{Nz_FeTop};
};

FeTop_6[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4006}; Layers{Nz_FeTop};
};

FeTop_7[] = Extrude {0., 0., Lz_FeTop} {
  Surface{4007}; Layers{Nz_FeTop};
};


Mesh 3;



// groups
Physical Surface("FSInterface") = {4005,4006,4007};

Physical Surface("FeTop_Down") = {4001,4002,4003,4004};
Physical Surface("FeTop_Up") = {FeTop_1[0],FeTop_2[0],FeTop_3[0],FeTop_4[0],FeTop_5[0],FeTop_6[0],FeTop_7[0]};
Physical Surface("FeTop_Front") = {FeTop_1[9],FeTop_2[9],FeTop_3[3]};
Physical Surface("FeTop_End") = {FeTop_1[3],FeTop_2[3],FeTop_4[5]};
Physical Surface("FeTop_Sides") = {FeTop_1[2],FeTop_2[2]};

Physical Volume("FeTop") = {FeTop_1[1],FeTop_2[1],FeTop_3[1],FeTop_4[1],FeTop_5[1],FeTop_6[1],FeTop_7[1]};
