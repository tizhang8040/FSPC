L = 110.0 ;
H = 3.0 ;

// Points ===============
Point(1) = {0,  0, 0};
Point(2) = {L,  0, 0};
Point(3) = {L,   H, 0};
Point(4) = {0  , H, 0};

// lines ================
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Curves ===============
Curve Loop(1) = {1, 2, 3, 4};

// Surfaces =============
Surface(1) = {1};

// Groups ===============
Physical Curve("FSInterface") = {3,4,2,1};
Physical Surface("Fluid") = {1};

// MESH REFINEMENT ======

//Transfinite Line {2, 4} = 6+1 Using Progression 1;
//Transfinite Line {1, 3} = 220+1 Using Progression 1;

//Transfinite Surface "*";

//Mesh.Algorithm = 9;