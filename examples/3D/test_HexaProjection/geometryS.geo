// parameters
Geometry.AutoCoherence = 0;	// keep duplicate elementary entities
// Fe-top
Lx_FeTop  = 120.;
Ly_FeTop  = 80.;
Lz_FeTop  = 1.5;
// Al-middle (solid)
Lx_Al = Lx_FeTop;
Ly_Al = Ly_FeTop;
Lz_Al = 3.;
// origin
X0 = 0.;
Y0 = 0.;
// Al-middle (solid)
Z0_Al = 0.0;
// Fe-top
Z0_FeTop = Z0_Al+Lz_Al;
//caracteristic lengths
// Al-middle (solid)
NxInt_Al = 50;
NxExt_Al = 30;
Nz_Al   = 2;
clInt_Al = Lx_Al/(NxInt_Al-1);
clExt_Al = Lx_Al/(NxExt_Al-1);
// Fe-top
NxInt_FeTop = 80;
NxExt_FeTop = 50;
Nz_FeTop   = 2;
clInt_FeTop = Lx_FeTop/(NxInt_FeTop-1);
clExt_FeTop = Lx_FeTop/(NxExt_FeTop-1);
//Printf("clInt_FeTop = %f", clInt_FeTop);



// Al-middle (solid)
// points
Point(1001) = {X0-Lx_Al/2, Y0-Ly_Al/2, Z0_Al, clExt_Al};
Point(1002) = {X0+Lx_Al/2, Y0-Ly_Al/2, Z0_Al, clExt_Al};
Point(1003) = {X0-Lx_Al/2, Y0,         Z0_Al, clInt_Al};
Point(1004) = {X0+Lx_Al/2, Y0,         Z0_Al, clInt_Al};
Point(1005) = {X0-Lx_Al/2, Y0+Ly_Al/2, Z0_Al, clExt_Al};
Point(1006) = {X0+Lx_Al/2, Y0+Ly_Al/2, Z0_Al, clExt_Al};

// lines
Line(1001) = {1001, 1002};
Line(1002) = {1002, 1004};
Line(1003) = {1004, 1003};
Line(1004) = {1003, 1001};
Line(1005) = {1004, 1006};
Line(1006) = {1006, 1005};
Line(1007) = {1005, 1003};

// curves
Curve Loop(1001) = { 1001, 1002, 1003, 1004};
Curve Loop(1002) = {-1003, 1005, 1006, 1007};

// surfaces
Plane Surface(1001) = {1001};
Plane Surface(1002) = {1002};

// 2D mesh
Recombine Surface{1001, 1002};

// volume
Al_1[] = Extrude {0., 0., Lz_Al} {
  Surface{1001}; Layers{Nz_Al}; Recombine;
};

Al_2[] = Extrude {0., 0., Lz_Al} {
  Surface{1002}; Layers{Nz_Al}; Recombine;
};


// Fe-top
// points
Point(2001) = {X0-Lx_FeTop/2, Y0-Ly_FeTop/2, Z0_FeTop, clExt_FeTop};
Point(2002) = {X0+Lx_FeTop/2, Y0-Ly_FeTop/2, Z0_FeTop, clExt_FeTop};
Point(2003) = {X0-Lx_FeTop/2, Y0,            Z0_FeTop, clInt_FeTop};
Point(2004) = {X0+Lx_FeTop/2, Y0,            Z0_FeTop, clInt_FeTop};
Point(2005) = {X0-Lx_FeTop/2, Y0+Ly_FeTop/2, Z0_FeTop, clExt_FeTop};
Point(2006) = {X0+Lx_FeTop/2, Y0+Ly_FeTop/2, Z0_FeTop, clExt_FeTop};

// lines
Line(2001) = {2001, 2002};
Line(2002) = {2002, 2004};
Line(2003) = {2004, 2003};
Line(2004) = {2003, 2001};
Line(2005) = {2004, 2006};
Line(2006) = {2006, 2005};
Line(2007) = {2005, 2003};

// curves
Curve Loop(2001) = { 2001, 2002, 2003, 2004};
Curve Loop(2002) = {-2003, 2005, 2006, 2007};

// surfaces
Plane Surface(2001) = {2001};
Plane Surface(2002) = {2002};

// 2D mesh
Recombine Surface{2001, 2002};

// volume
FeTop_1[] = Extrude {0., 0., Lz_FeTop} {
  Surface{2001}; Layers{Nz_FeTop}; Recombine;
};

FeTop_2[] = Extrude {0., 0., Lz_FeTop} {
  Surface{2002}; Layers{Nz_FeTop}; Recombine;
};


Mesh 3;


// Groups ===============
Physical Surface("Al_Down") = {1001,1002};
Physical Surface("Al_Up") = {Al_1[0],Al_2[0]};
Physical Surface("Al_Front") = {Al_1[5],Al_2[5]};
Physical Surface("Al_End") = {Al_1[3],Al_2[3]};
Physical Surface("Al_Sides") = {Al_1[2],Al_2[4]};

Physical Surface("FeTop_Down") = {2001,2002};
Physical Surface("FeTop_Up") = {FeTop_1[0],FeTop_2[0]};
Physical Surface("FeTop_Front") = {FeTop_1[5],FeTop_2[5]};
Physical Surface("FeTop_End") = {FeTop_1[3],FeTop_2[3]};
Physical Surface("FeTop_Sides") = {FeTop_1[2],FeTop_2[4]};

Physical Volume("Al") = {Al_1[1], Al_2[1]};
Physical Volume("FeTop") = {FeTop_1[1], FeTop_2[1]};

