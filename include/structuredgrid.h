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
    char* name;
    int np; //points in the poloidal magnetic line, elements number is np-1;
    int nr; //magnetic line number in radial direction, elements number is nr-1;

    // --- 2. Grid data ---
    CurveSet* zone_grid;

    // --- 3. Start & end tracing info ---
    double** start_points;   // [nr][2] starting point for tracing magnetic line
    double* guard_head;      // [nr]
    double* guard_end;       // [nr]
    double* pasmin;          // [nr]
    double* norm_pol_dist;  // distribution of points in the poloidal direction， from 0.0 to 1.0

    // --- 4. Boundary info ---
    Curve* first_boundary;
    Curve* second_boundary;

    // --- 5. Target info ---
    Curve* target_curve;
} Zone;

// Create a new Zone
Zone* allocate_zone();

// Load all geometry info from input file (returns 0 if success, 1 if fail)
int load_zone_from_file(Zone* z, const char* filename);

// Load individual curve components
int load_zone_first_boundary(Zone* z, Curve* first_boundary);
int load_zone_second_boundary(Zone* z, Curve* second_boundary);
int load_zone_target_curve(Zone* z, Curve* target_curve);

//free the zone
void free_zone(Zone** z);

#endif