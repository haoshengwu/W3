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
    int npol;  // Number of poloidal points (inner loop / fast index)
    int nrad;  // Number of radial points (outer loop / slow index)

    // Grid points stored in row-major order:
    //   - radial index ir = 0 to nrad-1 (outer loop)
    //   - poloidal index ip = 0 to npol-1 (inner loop)
    // Access: points[ir * npol + ip]
    GridPoint* points;
} TwoDimGrid;

TwoDimGrid* create_2Dgrid(int npol, int nrad);
void free_2Dgrid(TwoDimGrid* grid);

GridPoint get_point_2Dgrid(const TwoDimGrid* grid, int ip, int ir);
double get_x_2Dgrid(const TwoDimGrid* grid, int ip, int ir);
double get_y_2Dgrid(const TwoDimGrid* grid, int ip, int ir);
void set_point_2Dgrid(TwoDimGrid* grid, int ip, int ir, double x, double y);


//use CARRE algorithm to generate grid points
void generate_CARRE_2Dgrid(TwoDimGrid* grid,
                           GridZone* gridzone,
                           ode_function* func,
                           ode_solver* solver);

void generate_EMC3_2Dgrid(TwoDimGrid* grid,
                          GridZone* gridzone,
                          ode_function* func,
                          ode_solver* solver);



#endif