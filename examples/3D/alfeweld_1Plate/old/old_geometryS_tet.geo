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
hmax_steel              = 5.0;//3.0;
hmax_steel_tool         = 2.0;//1.0;
hmax_steel_FSI          = 1.25;//0.75;

hmax_backingPlate       = 4.5;//3.5;
hmax_backingPlate_FSI   = 2.0;//1.5;


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
             FINAL GROUPS
   ##################################### */

// all plates
Physical Surface("Top") = {6,8};
Physical Surface("FSInterface") = {7,1008,2005,2006,2007,2008};
Physical Volume("Solid") = {1};

