#ifndef TWODIMGRIDdGEN_H
#define TWODIMGIRDGEN_H

#include "gridzoneinfo.h"
#include "ode.h"
#include "equilibrium.h"
#include "sepdistribution.h"
#include "datastructure.h"


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  DLListNode* head;
  bool reverse;
} DLListWithOptions;

Curve* connect_DLList_for_curve(DLListWithOptions* list, int n); 


//GridZone is used to create the grid of TwoDimGrid!
typedef struct 
{
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
    // --- 3. End curve ---
    Curve* end_curve;

    // --- 4. first boundary ---
    bool first_bnd;
    Curve* first_bnd_curve;
    Curve* first_gridpoint_curve;

    // --- 5. Second boundary ---
    // used for multiple X-points situations. e.g. snowflakes
    bool sec_bnd;
    Curve* sec_bnd_curve;
    Curve* sec_gridpoint_curve;
} GridZone;


GridZone* create_sn_CARRE2D_GridZone(GridZoneInfo* gridzoneinfo, SepDistStr* sepdist);

void free_GridZone(GridZone* gridzone);

typedef struct {
    double x, y;  // Physical coordinates of a Grid point
} GridPoint;


//TwoDimGrid is written in one dimention to have fast access 
typedef struct {
    int npol;       // Logical number of poloidal points (columns)
    int nrad;       // Logical number of radial points (rows)
    int cap_npol;   // Allocated number of poloidal points (row stride)
    int cap_nrad;   // Allocated number of radial points (column count)
    int offset_rad; // Row offset from physical memory origin
    int offset_pol; // Column offset from physical memory origin

    // Grid points stored in row-major order:
    //   - ir = 0 to nrad-1 maps to offset_rad + ir
    //   - ip = 0 to npol-1 maps to offset_pol + ip
    // Access: points[(offset_rad + ir) * cap_npol + (offset_pol + ip)]
    GridPoint* points;
} TwoDimGrid;

TwoDimGrid* create_2Dgrid_default(int npol, int nrad);
void free_2Dgrid(TwoDimGrid* grid);

GridPoint* get_point_2Dgrid(TwoDimGrid* g, int ir, int ip);
double get_x_2Dgrid(const TwoDimGrid* g, int ir, int ip);
double get_y_2Dgrid(const TwoDimGrid* g, int ir, int ip);
void set_point_2Dgrid(TwoDimGrid* g, int ir, int ip, double x, double y);


//use CARRE algorithm to generate grid points for on zone
// TwoDimGrid* grid is the grid
// GridZone* grizone store the necessary data for grid generation
// ode_function* func and ode_solver* solver are used for line tracing
void generate_CARRE_2Dgrid_default(TwoDimGrid* grid,
                                   GridZone* grizone,
                                   ode_function* func,
                                   ode_solver* solver);

void generate_EMC3_2Dgrid(TwoDimGrid* grid,
                          GridZoneInfo* gridzoneinfo,
                          ode_function* func,
                          ode_solver* solver);



#endif