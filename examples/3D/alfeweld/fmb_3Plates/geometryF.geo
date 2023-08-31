// parameters
Geometry.AutoCoherence = 0;	// keep duplicate elementary entities
// geometry in [mm]
// Al-middle (fluid)
Lx_AlF = 110.;
Ly_AlF = 20.;
Lz_AlF = 3.;
gx_AlF = 12.5;			// gap in x direction between left edges of AlS-AlF
// origin
X0 = 0.;
Y0 = 0.;
Z0_AlF = 5.;
//caracteristic lengths
// Al-middle (fluid)
cl_AlF = 1.2;



// Al-middle (fluid)
// points
Point(3001) = {X0-Lx_AlF/2, Y0-Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3002) = {X0+Lx_AlF/2, Y0-Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3003) = {X0+Lx_AlF/2, Y0+Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3004) = {X0-Lx_AlF/2, Y0+Ly_AlF/2, Z0_AlF,        cl_AlF};
Point(3005) = {X0-Lx_AlF/2, Y0-Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};
Point(3006) = {X0+Lx_AlF/2, Y0-Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};
Point(3007) = {X0+Lx_AlF/2, Y0+Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};
Point(3008) = {X0-Lx_AlF/2, Y0+Ly_AlF/2, Z0_AlF+Lz_AlF, cl_AlF};

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


Mesh 3;


// Groups ===============
Physical Surface("FSInterface") = {3001, 3002, 3003, 3004, 3005, 3006};
Physical Volume("AlF") = {3001};

