#ifndef TWODIMGRIDdGEN_H
#define TWODIMGIRDGEN_H

#include "gridzone.h"
#include "ode.h"
#include "equilibrium.h"


#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double x, y;  // Physical coordinates of a Grid point
} GridPoint;

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


//use CARRE algorithm to generate grid points
void generate_CARRE_2Dgrid(TwoDimGrid* grid,
                           GridZoneInfo* gridzoneinfo,
                           ode_function* func,
                           ode_solver* solver);

void generate_EMC3_2Dgrid(TwoDimGrid* grid,
                          GridZoneInfo* gridzoneinfo,
                          ode_function* func,
                          ode_solver* solver);



#endif