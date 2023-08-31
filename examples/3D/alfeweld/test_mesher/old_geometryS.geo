L = 20 ;
H = 3 ;
Hs= 1.5 ; 
W = 10;

// Points ===============
Point(1) = {0,  H , 0};
Point(2) = {L,  H , 0};
Point(3) = {L,  H + Hs, 0};
Point(4) = {0  ,H + Hs, 0};

Point(5) = {0,  H ,   W};
Point(6) = {L,  H ,   W};
Point(7) = {L,  H + Hs, W};
Point(8) = {0  ,H + Hs, W};

// lines ================
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {6, 5};
Line(6) = {5, 8};
Line(7) = {8, 7};
Line(8) = {7, 6};

Line(9)  = {6, 2};
Line(10) = {3, 7};
Line(11) = {8, 4};
Line(12) = {1, 5};


// Curves ===============
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {-8, -10, -2, -9};
Curve Loop(4) = {-1, 12, -5, 9};
Curve Loop(5) = {-4, -11, -6, -12};
Curve Loop(6) = {-7, 11, -3, 10};

// Surfaces =============
Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};
Surface(4) = {4};
Surface(5) = {5};
Surface(6) = {6};

// volume ===============
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Groups ===============
Physical Surface("Top") = {6};
Physical Surface("Clamped") = {1,2,3,5};
Physical Surface("FSInterface") = {4};
Physical Volume("Solid") = {1};