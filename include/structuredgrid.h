#ifndef STRUCTUREDGRID_H
#define STRUCTUREDGRID_H
#include <stdio.h>
#include "curve.h"

typedef struct {
    // --- 1. Basic information ---
    char* name;
    int nr; //magnetic line number in radial direction, elements number is nr-1;
    // --- 2. Start & end tracing info ---
    double* start_point_R;  // [nr]starting point R coordinate for tracing magnetic line
    double* start_point_Z;  // [nr]starting point R coordinate for tracing magnetic line
    double* guard_start;    // [nr]
    double* guard_end;      // [nr]
    double* pasmin;         // [nr]

    //first line distribution
    int npoint; //points in the poloidal magnetic line, elements number is np-1;
    double* norm_pol_dist;  // [np]distribution of points in the poloidal direction， from 0.0 to 1.0

    // --- 3. Boundary info ---
    int n_boundary;
    Curve* first_boundary;
    //TODO
    //Curve* second_boundary;
    Curve* end_curve;

    // --- 5. Grid data ---
    CurveSet* grid_curveset;
} GridZone;

// Create a new GridZone
GridZone* allocate_GridZone();


//free the GridZone
void free_GridZone(GridZone** z);

//Write the GridZone to input file for mesh generation;
void write_input_from_GridZone(GridZone* gridzone, const char* filename);

//create the GridZone from a input file for mesh generation;
GridZone* load_GridZone_from_input(const char* filename);

//print GirdZone info for test. first_boundary and curveset do not print;
void print_GridZone(GridZone* gridzone);

#endif