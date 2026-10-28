// =============================================================
// Cylinder flow — FIXED physical group IDs (no ID conflicts)
// Domain: [-15, 25] x [-8, 8], D=1 cylinder at origin
// Re=100 (u_inf=1, D=1, nu=0.01)
// =============================================================

lc_cyl    = 0.15;
lc_wake   = 0.5;
lc_far    = 1.5;

// Domain corners
Point(1)  = {-15, -8, 0, lc_far};
Point(2)  = { 25, -8, 0, lc_far};
Point(3)  = { 25,  8, 0, lc_far};
Point(4)  = {-15,  8, 0, lc_far};

// Cylinder (radius=0.5)
Point(5)  = { 0,    0,   0, lc_cyl};
Point(6)  = { 0.5,  0,   0, lc_cyl};
Point(7)  = { 0,    0.5, 0, lc_cyl};
Point(8)  = {-0.5,  0,   0, lc_cyl};
Point(9)  = { 0,   -0.5, 0, lc_cyl};

// Outer boundary lines
Line(1) = {1, 2};   // bottom
Line(2) = {2, 3};   // right (outlet)
Line(3) = {3, 4};   // top
Line(4) = {4, 1};   // left (inlet)

// Cylinder arcs
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};

// =============================================================
// PHYSICAL GROUPS — all unique IDs, surface gets highest ID
// Curve tags: 1,2,3,4  (1D boundaries)
// Surface tag: 100      (2D domain — use a high number, never 1-4)
// =============================================================
Physical Curve("inlet",    1) = {4};
Physical Curve("outlet",   2) = {2};
Physical Curve("walls",    3) = {1, 3};
Physical Curve("cylinder", 4) = {5, 6, 7, 8};
Physical Surface("fluid", 100) = {1};   // <-- tag 100, NOT 1

// Refinement fields
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling   = 100;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_cyl;
Field[2].SizeMax  = lc_far;
Field[2].DistMin  = 0.0;
Field[2].DistMax  = 6.0;

Field[3] = Box;
Field[3].XMin = -0.5;  Field[3].XMax = 25.0;
Field[3].YMin = -2.5;  Field[3].YMax =  2.5;
Field[3].VIn  = lc_wake;
Field[3].VOut = lc_far;

Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;

Mesh.Algorithm    = 5;
Mesh.RecombineAll = 0;
Mesh.Smoothing    = 5;
Mesh.ElementOrder = 1;
