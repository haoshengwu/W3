#ifndef STRUCTUREDGRID_H
#define STRUCTUREDGRID_H
#include <stdio.h>
#include "curve.h"

typedef struct {
    // --- 1. Basic information ---
    char* name;
    int np; //points in the poloidal magnetic line, elements number is np-1;
    int nr; //magnetic line number in radial direction, elements number is nr-1;

    // --- 2. Grid data ---
    CurveSet* grid_curveset;

    // --- 3. Start & end tracing info ---
    double* start_point_R;  // [nr]starting point R coordinate for tracing magnetic line
    double* start_point_Z;  // [nr]starting point R coordinate for tracing magnetic line
    double* guard_start;    // [nr]
    double* guard_end;      // [nr]
    double* pasmin;         // [nr]
    double* norm_pol_dist;  // distribution of points in the poloidal direction， from 0.0 to 1.0

    // --- 4. Boundary info ---
    int n_boundary;
    Curve* first_boundary;
    Curve* second_boundary;

    // --- 5. Target info ---
    Curve* target_curve;
} GridZone;

// Create a new GridZone
GridZone* allocate_GridZone();

// Load all geometry info from input file (returns 0 if success, 1 if fail)
int load_GridZone_from_file(GridZone* z, const char* filename);

// Load individual curve components
int load_GridZone_first_boundary(GridZone* z, Curve* first_boundary);
int load_GridZone_second_boundary(GridZone* z, Curve* second_boundary);
int load_GridZone_target_curve(GridZone* z, Curve* target_curve);

//free the GridZone
void free_GridZone(GridZone** z);

#endif