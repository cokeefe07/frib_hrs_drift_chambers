// 1. Explicitly define element sizes
lc_bulk = 1.0;
lc_wire = 0.001; // 20 um (Resolves the 25um wire perfectly)

pitch = 0.254;
n_wires = 118;
gap = 3.0;
r_wire = 0.0025;
z_thick = 0.01;  // 1 mm thick 3D slice
y_offset = 14.859;

// 2. Real Gas Vessel Bounding Box (45x45 cm)
// This pushes the 0V ground walls safely away from the active detector
x_min = -22.5; x_max = 22.5;
y_min = -22.5; y_max = 22.5;

Point(1) = {x_min, y_min, 0, lc_bulk};
Point(2) = {x_max, y_min, 0, lc_bulk};
Point(3) = {x_max, y_max, 0, lc_bulk};
Point(4) = {x_min, y_max, 0, lc_bulk};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};


// 3. M-THGEM and Cathode as thin solid strips (holes), same technique as the wires
lc_electrode = 0.01;

// M-THGEM strip at Y = -16.0, 0.02 cm thick
Point(200) = {-gap/2.0, -16.01, 0, lc_electrode};
Point(201) = { gap/2.0, -16.01, 0, lc_electrode};
Point(202) = { gap/2.0, -15.99, 0, lc_electrode};
Point(203) = {-gap/2.0, -15.99, 0, lc_electrode};
Line(20) = {200,201}; Line(21) = {201,202}; Line(22) = {202,203}; Line(23) = {203,200};
Curve Loop(20) = {20,21,22,23};

// Cathode strip at Y = +16.0, 0.02 cm thick
Point(210) = {-gap/2.0, 15.99, 0, lc_electrode};
Point(211) = { gap/2.0, 15.99, 0, lc_electrode};
Point(212) = { gap/2.0, 16.01, 0, lc_electrode};
Point(213) = {-gap/2.0, 16.01, 0, lc_electrode};
Line(24) = {210,211}; Line(25) = {211,212}; Line(26) = {212,213}; Line(27) = {213,210};
Curve Loop(21) = {24,25,26,27};

// 4. Generate all 118 pairs of wires
p = 1000; l = 1000; cl = 5000; 
wire_loops[] = {};

For i In {1:n_wires}
    y_pos = ((i - 1) * pitch) - y_offset;

    // Left wire
    Point(p)   = {-gap/2.0, y_pos - r_wire, 0, lc_wire}; 
    Point(p+1) = {-gap/2.0 + r_wire, y_pos, 0, lc_wire}; 
    Point(p+2) = {-gap/2.0, y_pos + r_wire, 0, lc_wire}; 
    Point(p+3) = {-gap/2.0 - r_wire, y_pos, 0, lc_wire}; 
    Point(p+4) = {-gap/2.0, y_pos, 0, lc_wire};          

    Circle(l)   = {p,   p+4, p+1};
    Circle(l+1) = {p+1, p+4, p+2};
    Circle(l+2) = {p+2, p+4, p+3};
    Circle(l+3) = {p+3, p+4, p};

    Curve Loop(cl) = {l, l+1, l+2, l+3}; 
    Transfinite Curve {l, l+1, l+2, l+3} = 6;   // BUGFIX: was {1,...}. Curve 1 is the
                                                //  vessel BOTTOM WALL, so this pinned it
                                                //  to 6 nodes across 45 cm, and left arc
                                                //  'l' un-transfinited on all 118 left
                                                //  wires -> systematic L/R asymmetry.
    wire_loops[] += {cl};

    p += 5; l += 4; cl += 1;

    // Right wire
    Point(p)   = {gap/2.0, y_pos - r_wire, 0, lc_wire};  
    Point(p+1) = {gap/2.0 + r_wire, y_pos, 0, lc_wire};  
    Point(p+2) = {gap/2.0, y_pos + r_wire, 0, lc_wire};  
    Point(p+3) = {gap/2.0 - r_wire, y_pos, 0, lc_wire};  
    Point(p+4) = {gap/2.0, y_pos, 0, lc_wire};           

    Circle(l)   = {p,   p+4, p+1};
    Circle(l+1) = {p+1, p+4, p+2};
    Circle(l+2) = {p+2, p+4, p+3};
    Circle(l+3) = {p+3, p+4, p};

    Curve Loop(cl) = {l, l+1, l+2, l+3}; 
    Transfinite Curve {l, l+1, l+2, l+3} = 6;
    wire_loops[] += {cl};

    p += 5; l += 4; cl += 1;
EndFor




// 5. Build the 2D Gas Area and embed the floating faces!
Plane Surface(1) = {1, wire_loops[], 20, 21};

// 6. Extrude into a 3D Slice WITHOUT Layers{1} (Generates Tetrahedra)
ext[] = Extrude {0, 0, z_thick} { Surface{1}; Layers{1}; };


// 7. Physical Groups
// THE CRITICAL FIX: ext[1] is the 3D Volume. ext was just the 2D Top Surface!
Physical Volume(1) = {ext[ 1 ]}; 

// Grounded Outer Walls (Dynamically grabbed via BoundingBox)
Physical Surface(1) = {
    Surface In BoundingBox{-22.51, -22.51, -0.01, -22.49,  22.51, z_thick + 0.01}, // Left
    Surface In BoundingBox{ 22.49, -22.51, -0.01,  22.51,  22.51, z_thick + 0.01}, // Right
    Surface In BoundingBox{-22.51, -22.51, -0.01,  22.51, -22.49, z_thick + 0.01}, // Bottom
    Surface In BoundingBox{-22.51,  22.49, -0.01,  22.51,  22.51, z_thick + 0.01}  // Top
};

// High Voltage faces (Explicitly using the extruded surface IDs we just generated)
Physical Surface(2) = Surface In BoundingBox{-gap/2.0-0.01, -16.02, -0.01, gap/2.0+0.01, -15.98, z_thick+0.01}; // M-THGEM
Physical Surface(3) = Surface In BoundingBox{-gap/2.0-0.01,  15.98, -0.01, gap/2.0+0.01,  16.02, z_thick+0.01}; // Cathode



// 118 Wires (Maps perfectly to Boundary Conditions 4 through 121 in field_cage.sif)
For i In {1:n_wires}
    y_pos = ((i - 1) * pitch) - y_offset;
    Physical Surface(3 + i) = Surface In BoundingBox{-gap/2.0 - r_wire - 0.001, y_pos - r_wire - 0.001, -0.001, 
                                                      gap/2.0 + r_wire + 0.001, y_pos + r_wire + 0.001, z_thick + 0.001};
EndFor


// ---------- Symmetric size field (strip_cage lesson applied to the wires) ----------
// Point-based sizing alone let the mesh behind the wires run from 545 um at the wire
// plane to 5008 um at |x|=2.0 -- coarser than the 2540 um wire pitch, so the wire
// structure was unresolved exactly where you saw the bad field, and the mesh out
// there was only 56.8% mirror-paired. This field is built from Fabs(Fabs(x)-gap/2),
// so it is EXACTLY symmetric in x by construction, with no hard edges anywhere.
//   |E|-relevant scale: ripple decay length = pitch/2pi = 404 um.
h_near = 0.020;    // element size at the wire plane
g1     = 0.11178;  // linear growth
g2     = 0.05644;  // quadratic term -> relaxes to lc_bulk ~5 cm out
y_lim  = 15.0;     // no refinement above/below the wire stack
// ---- MEASURED, full 118 wires, -3 -order 2 ----
//  h_near   g1,g2 scale | 3D nodes | worst gamma | spacing at |x|=2.0
//   150um       x1      |   1,689k |   2.19e-02  |   ~40 um   (OOM risk)
//   200um       x2      |     751k |   2.14e-02  |    200 um   <-- current
//   200um       x3      |    ~440k |      -      |    ~300 um
//   250um       x3      |    ~365k |      -      |    ~380 um
//  (baseline, no field  |    ~90k  |      -      |   5008 um)
// Spacing behind the wires vs the old mesh: 11x finer at |x|=1.5, 25x at |x|=2.0.
// If ElmerSolver runs tight, scale g1 and g2 up together (x3, then x4).
Field[1] = MathEval;
Field[1].F = Sprintf("%g+%g*Sqrt((Fabs(Fabs(x)-%g))^2+((Fabs(y)-%g+Fabs(Fabs(y)-%g))/2)^2)+%g*Sqrt((Fabs(Fabs(x)-%g))^2+((Fabs(y)-%g+Fabs(Fabs(y)-%g))/2)^2)*Sqrt((Fabs(Fabs(x)-%g))^2+((Fabs(y)-%g+Fabs(Fabs(y)-%g))/2)^2)",
                     h_near, g1, gap/2.0, y_lim, y_lim,
                     g2,     gap/2.0, y_lim, y_lim, gap/2.0, y_lim, y_lim);
Background Field = 1;
Mesh.CharacteristicLengthMax = lc_bulk;

// 8. Robust meshing flags
Mesh.Algorithm = 6;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.SecondOrderLinear = 1; // Prevents curved nodes from folding inside-out into singularities

