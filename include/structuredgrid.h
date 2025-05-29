#ifndef STRUCTUREDGRID_H
#define STRUCTUREDGRID_H
#include <stdio.h>


// Represents a single 2D curve with n_point points
typedef struct {
    size_t n_point;         // Number of (x, y) points in the curve
    double** points;      // Array of points [n_point][2], each point is [x, y]
} Curve;

// Represents a set of curves
typedef struct {
    size_t n_curve;      // Number of curves in the set
    Curve** curves;       // Array of Curve structures
} CurveSet;

// Allocate memory for a single curve with 'n_point' (x,y) pairs
Curve* create_curve(size_t n_point);

// Free memory of a single curve
void free_curve(Curve* c);

// Allocate memory for a curve set with 'n_curve' curves,
// each containing 'n_point' points
CurveSet* create_curveset(size_t n_curve, size_t n_point);

// Free memory of the entire curve set (all curves and their points)
void free_curveset(CurveSet* cs);

typedef struct {
    // --- 1. Basic information ---
    char* name;              // Zone name (e.g., PFR123, SOL123, CORE)
    int np;                  // Number of magnetic surfaces (radial direction)
    int nr;                  // Number of points along each surface (poloidal)

    // --- 2. Grid data ---
    // zone_grid is a set of cureve which is consist of [nr] curves,
    // and each curve is consist of [np] points.
    CurveSet* zone_grid;

    // --- 3. Start & end tracing info ---
    double** start_points;   // [nr][2] Starting points for field line tracing
    double* guard_head;      // [nr] Guard length at the start of each line
    double* guard_end;       // [nr] Guard length at the end of each line
    double* pasmin;          // [nr] Minimum parallel arc length (pasmin) for each line
    double distribution[2];  // deltp1 and deltpn from CARRE

    // --- 4. Boundary info ---
    int n_inner_boundary;    
    double** inner_boundary;  // [n_inner_boundary][2]

    int n_outer_boundary;    
    double** outer_boundary;  // [n_outer_boundary][2]

    // --- 5. Target info ---
    int n_target;           
    double** target_points;   // [n_target][2]
    // --- 6. Miscellaneous ---
} Zone;








#endif