L = 40 ;
H = 3 ;
Hs= 1.5 ; 
W = 30;

HB= 10;

LL= 20;
LR= 40;
WR= 20;
WL= 20;

/* #####################################
            ELEMENT SIZE CONTROL
   ##################################### */
hmax_steel              = 10.0;//3.0;
hmax_steel_tool         = 6.0;//1.0;
hmax_steel_FSI          = 4.0;//0.75;

hmax_backingPlate       = 10.0;//3.5;
hmax_backingPlate_FSI   = 5.5;//1.5;


/* #####################################
               STEEL PLATE
   ##################################### */

// Points Steel OUT ==================
Point(1) = {-LL ,  H     , -WL , hmax_steel};
Point(2) = {L+LR,  H     , -WL , hmax_steel};
Point(3) = {L+LR,  H + Hs, -WL , hmax_steel}; 
Point(4) = {-LL ,  H + Hs, -WL , hmax_steel};
Point(5) = {-LL ,  H     , WR+W, hmax_steel};
Point(6) = {L+LR,  H     , WR+W, hmax_steel};
Point(7) = {L+LR,  H + Hs, WR+W, hmax_steel};
Point(8) = {-LL ,  H + Hs, WR+W, hmax_steel};

// Points Steel IN  ==================
Point(9)  = {0 ,   H      , 0, hmax_steel_FSI};
Point(10) = {0 ,   H      , W, hmax_steel_FSI};
Point(11) = {L ,   H      , W, hmax_steel_FSI};
Point(12) = {L ,   H      , 0, hmax_steel_FSI};
Point(13) = {0 ,   H + Hs , 0, hmax_steel_tool};
Point(14) = {0 ,   H + Hs , W, hmax_steel_tool};
Point(15) = {L ,   H + Hs , W, hmax_steel_tool};
Point(16) = {L ,   H + Hs , 0, hmax_steel_tool};

// lines Steel OUT ===================
Line(1) = {1, 2};  Line(2) = {2, 3};
Line(3) = {3, 4};  Line(4) = {4, 1};
Line(5) = {6, 5};  Line(6) = {5, 8};
Line(7) = {8, 7};  Line(8) = {7, 6};
Line(9)  = {6, 2}; Line(10) = {3, 7};
Line(11) = {8, 4}; Line(12) = {1, 5};

// lines Steel IN ====================
Line(13) = {10, 9};  Line(14) = {9, 12};
Line(15) = {12, 11}; Line(16) = {11, 10};
Line(17) = {13, 14}; Line(18) = {14, 15};
Line(19) = {15, 16}; Line(20) = {16, 13};

// Curves ============================
Curve Loop(1) = {-1, -2, -3, -4};
Curve Loop(2) = {-5, -6, -7, -8};
Curve Loop(3) = {8, 10, 2, 9};
Curve Loop(4) = {1, -12, 5, -9};
Curve Loop(5) = {4, 11, 6, 12};
Curve Loop(6) = {7, -11, 3, -10};
Curve Loop(7) = {13, 14, 15, 16};
Curve Loop(8) = {17, 20, 19, 18};

// Surfaces ==========================
Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};
Plane Surface(4) = {4,7};
Surface(5) = {5};
Plane Surface(6) = {6,8};
Surface(7) = {7};
Surface(8) = {8};

// volume ===============
Surface Loop(1) = {1, 2, 3, 4, 7, 5, 6, 8};
Volume(1) = {1};

/* #####################################
               BACKING PLATE
   ##################################### */

// Points Steel OUT ==================
Point(1001) = {-LL ,  -HB     , -WL , hmax_backingPlate};
Point(1002) = {L+LR,  -HB     , -WL , hmax_backingPlate};
Point(1003) = {L+LR,  -HB + HB, -WL , hmax_backingPlate}; 
Point(1004) = {-LL ,  -HB + HB, -WL , hmax_backingPlate};
Point(1005) = {-LL ,  -HB     , WR+W, hmax_backingPlate};
Point(1006) = {L+LR,  -HB     , WR+W, hmax_backingPlate};
Point(1007) = {L+LR,  -HB + HB, WR+W, hmax_backingPlate};
Point(1008) = {-LL ,  -HB + HB, WR+W, hmax_backingPlate};


// Points Steel IN  ==================
Point(1013) = {0 ,   -HB + HB , 0, hmax_backingPlate_FSI};
Point(1014) = {0 ,   -HB + HB , W, hmax_backingPlate_FSI};
Point(1015) = {L ,   -HB + HB , W, hmax_backingPlate_FSI};
Point(1016) = {L ,   -HB + HB , 0, hmax_backingPlate_FSI};

// lines Steel OUT ===================
Line(1001) = {1001, 1002};  Line(1002) = {1002, 1003};
Line(1003) = {1003, 1004};  Line(1004) = {1004, 1001};
Line(1005) = {1006, 1005};  Line(1006) = {1005, 1008};
Line(1007) = {1008, 1007};  Line(1008) = {1007, 1006};
Line(1009) = {1006, 1002};  Line(1010) = {1003, 1007};
Line(1011) = {1008, 1004};  Line(1012) = {1001, 1005};

// lines Steel IN ====================
Line(1017) = {1013, 1014}; Line(1018) = {1014, 1015};
Line(1019) = {1015, 1016}; Line(1020) = {1016, 1013};


// Curves ============================
Curve Loop(1001) = {-1001, -1002, -1003, -1004};
Curve Loop(1002) = {-1005, -1006, -1007, -1008};
Curve Loop(1003) = {1008, 1010, 1002, 1009};
Curve Loop(1004) = {1001, -1012, 1005, -1009};
Curve Loop(1005) = {1004, 1011, 1006, 1012};
Curve Loop(1006) = {1007, -1011, 1003, -1010};
Curve Loop(1008) = {1017, 1020, 1019, 1018};


// Surfaces ==========================
Surface(1001) = {1001};
Surface(1002) = {1002};
Surface(1003) = {1003};
Surface(1004) = {1004};
Surface(1005) = {1005};
Plane Surface(1006) = {1006,1008};
Surface(1008) = {1008};

// volume ===============
Surface Loop(2) = {1001, 1002, 1003, 1004, 1005, 1006, 1008};
Volume(2) = {2};   



/* #####################################
             aluminium PLATE
   ##################################### */

Line(2001) = {1, 1004};  
Line(2002) = {2, 1003};  
Line(2003) = {5, 1008};  
Line(2004) = {6, 1007};  
Line(2005) = {9, 1013};  
Line(2006) = {12,1016};
Line(2007) = {10,1014};
Line(2008) = {11,1015};  

Curve Loop(2001) = {-12, -2003, -1011, 2001};
Curve Loop(2002) = {5 , 2003, 1007, -2004};
Curve Loop(2003) = {-9, 2004, -1010,-2002};
Curve Loop(2004) = {1 , 2002, 1003, -2001};
Curve Loop(2005) = {-13,2007,-1017, -2005};
Curve Loop(2006) = {-16,2008,-1018, -2007};
Curve Loop(2007) = {-15,2006,-1019, -2008};
Curve Loop(2008) = {-14,2005,-1020, -2006};

Surface(2001) = {2001};
Surface(2002) = {2002};
Surface(2003) = {2003};
Surface(2004) = {2004};
Surface(2005) = {2005};
Surface(2006) = {2006};
Surface(2007) = {2007};
Surface(2008) = {2008};

// inverted surfaces
Curve Loop(3004) = {-1, 12, -5, 9};
Curve Loop(3007) = {-13, -14, -15, -16};
Plane Surface(3004) = {3004,3007};

Curve Loop(3006) = {-1007, 1011, -1003, 1010};
Curve Loop(3008) = {-1017, -1020, -1019, -1018};
Plane Surface(3006) = {3006,3008};

Surface Loop(3) = {3004,2001,2002,2003,2004,2005,2006,2007,2008,3006};

Volume(3) = {3};

Mesh.SubdivisionAlgorithm = 2;


/* #####################################
             FINAL GROUPS
   ##################################### */

// all plates
Physical Surface("Top") = {8};
Physical Surface("Clamped") = {1,2,3,5,6,2001,2002,2003,2004,1001,1002,1003,1005};
Physical Surface("FSInterface") = {7,1008,2005,2006,2007,2008};
Physical Surface("Bottom") = {1004};
Physical Volume("Solid") = {1,2,3};

