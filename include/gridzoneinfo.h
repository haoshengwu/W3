#ifndef GRIDZONE_H
#define GRIDZONE_H
#include "curve.h"


//Information of GridZone, which is mainly from DivGeo *.trg fiel.
typedef struct {
    // --- 1. Basic information ---
    char* topo;
    char* name;
    // --- 2. Start & end tracing info ---
    int nr; //magnetic line number in radial direction, elements number is nr-1;
    double* start_point_R;  // [nr]starting point R coordinate for tracing magnetic line
    double* start_point_Z;  // [nr]starting point R coordinate for tracing magnetic line
    double* guard_start;    // [nr]
    double* guard_end;      // [nr]
    double* pasmin;         // [nr]

    // --- 3. Start and End curves ---
    Curve* start_curve;
    Curve* end_curve;

    // --- 4. first poloidal boundary info
    int n_polsegm1;
    int* xptidx1; //indicate index of 
    int* seplineidx1;
    int* segmidx1;
    //indicate whether segments need to be reserve or not when create the start curve
    // 1 means reverse, 0 means not;
    int* reverse_segm1;

    // --- 5. second poloidal boundary info
    int n_polsegm2;
    int* xptidx2; //indicate index of 
    int* seplineidx2;
    int* segmidx2;
    int* reverse_segm2;
} GridZoneInfo;

// Create a new GridZoneInfo
GridZoneInfo* allocate_GridZoneInfo();

//free the GridZoneInfo
void free_GridZoneInfo(GridZoneInfo** z);

//Write the GridZoneInfo to input file for mesh generation;
void write_GridZoneInfo(GridZoneInfo* gridzoneinfo, const char* filename);

//create the GridZoneInfo from a input file for mesh generation;
GridZoneInfo* load_GridZoneInfo_from_input(const char* filename);

//print GirdZone info for test. first_pol_points and curveset do not print;
void print_GridZoneInfo(GridZoneInfo* gridzoneinfo);

#endif