d = 0.01;
L = 110.0 ;
H = 3.0 ;
Hs= 1.5 ; 
Hb= 10  ;
dL= 12.5;

// Points ===============
Point(1) = {L, H, 0, d};
Point(2) = {0, H, 0, d};
Point(3) = {L, H+Hs, 0, d};
Point(4) = {0, H+Hs, 0, d};

Point(5) = {-dL , H+Hs , 0, d};
Point(6) = {-dL , H , 0, d};
Point(7) = {L+dL, H+Hs , 0, d};
Point(8) = {L+dL, H , 0, d};

Point(9)  = {-dL , 0 , 0, d};
Point(10) = { 0  , 0 , 0, d};
Point(11) = {L   , 0 , 0, d};
Point(12) = {L+dL, 0 , 0, d};

Point(13) = {-dL , -Hb , 0, d};
Point(14) = { 0  , -Hb , 0, d};
Point(15) = {L   , -Hb , 0, d};
Point(16) = {L+dL, -Hb , 0, d};

// lines ================
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Line(5)  = {5, 6};
Line(6)  = {6, 2};
Line(7)  = {4, 5};
Line(8)  = {1, 8};
Line(9)  = {8, 7};
Line(10) = {7, 3};

Line(11) = {6, 9};
Line(12) = {9, 10};
Line(13) = {10, 2};
Line(14) = {1, 11};
Line(15) = {11, 12};
Line(16) = {12, 8};

Line(17) = {9,  13};
Line(18) = {13, 14};
Line(19) = {14, 10};
Line(20) = {11, 15};
Line(21) = {15, 16};
Line(22) = {16, 12};

Line(23) = {10, 11};
Line(24) = {14, 15};

// Curves ===============
Curve Loop(1) = {-1, -2, -3, -4};
Curve Loop(2) = {5, 6, 2, 7};
Curve Loop(3) = {8, 9, 10, 4};
Curve Loop(4) = {11, 12, 13, -6};
Curve Loop(5) = {14, 15, 16, -8};
Curve Loop(6) = {17, 18, 19, -12};
Curve Loop(7) = {20, 21, 22, -15};
Curve Loop(8) = {-19, 24, -20, -23};

// Surfaces =============
Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};
Surface(4) = {4};
Surface(5) = {5};
Surface(6) = {6};
Surface(7) = {7};
Surface(8) = {8};

// Groups ===============
Physical Curve("Top") = {3,7,10};
Physical Curve("Clamped") = {5,11,17,9,16,22};
Physical Curve("Bottom") = {18,24,21};
Physical Curve("FSInterface") = {1,13,14,23};
Physical Surface("Solid_1") = {1,2,3,4,5,6,7,8};

// MESH 
Transfinite Line {1, 3} = 400 Using Progression 1;
Transfinite Line {2, 4, 5, 9} = 10 Using Progression 1;

Transfinite Line {7, 6, 12, 18, 8, 10, 15, 21} = 10 Using Progression 1;
Transfinite Line {11, 13, 14, 16} = 20 Using Progression 1;

Transfinite Line {17, 19, 20, 22} = 15 Using Progression 1;

Transfinite Line {23,24} = 100 Using Progression 1;


Recombine Surface "*";
Transfinite Surface "*";

Mesh.Binary=1;
