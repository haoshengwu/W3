#ifndef TFI2DGRIDGEN_H
#define TFI2DGRIDGEN_H

#include "curve.h"
#include "twodimgridgen.h"
#include "config.h"


/*
Generate the 2d strucutred grid by the Transfinite Interpolation (TFI) 
grid is the input of a two initialized grid
cur_l is the grid point curve of left boundary
cur_r is the grid point curve of right boundary
cur_b is the grid point curve of bottom boundary
cur_t is the grid point curve of bottom boundary
distb_cur_l is the normalized distribution of grid point in cur_l, nleft is the size
distb_cur_r is the normalized distribution of grid point in cur_r, nright is the size
distb_cur_b is the normalized distribution of grid point in cur_b, nbottom is the size
distb_cur_t is the normalized distribution of grid point in cur_t, top is the size

Following is the schematic view
                    cur_t ->
        :-----------------------------:
      ^ :                             :   ^
      | :                             :   |
   cur_l:                             : cur_r
        :                             :
        :                             :
        :-----------------------------:
                    cur_b ->
*/

void check_tfi_boundary(Curve* cur_l, Curve* cur_r, Curve* cur_b, Curve* cur_t);

void generate_2Dgrid_default_TFI(TwoDimGrid* grid,
                                 Curve* cur_b, Curve* cur_t, Curve* cur_l, Curve* cur_r,
                                 const double* distb_cur_b, int nbottom,
                                 const double* distb_cur_t, int ntop,
                                 const double* distb_cur_l, int nleft, 
                                 const double* distb_cur_r, int nright);


//fixed the grid points which from the 2d grid from TFI.
//in some situation, the grid points is out ouf boundary.
//The funcion can fixed the points.
void optimized_neu_2Dgrid(TwoDimGrid* grid);

int try_neighbor_averaging(TwoDimGrid* grid, Curve** boundaris, int n_boundary, int max_iter, double drmin);


//if the distance between [i][j-1] [i][j] < drmin then average the point according the neighbor
void try_avoid_drmin(TwoDimGrid* grid, double drmin);

//Count how many points in the grid are inside the boundary. 
//the points are [1:end-1][1:end-1].
int count_points_inside(TwoDimGrid* grid, Curve** boundary, int n_boundary);


//use LaplaceSmoothing to smooth the 2d grid return the iteration number
int LaplaceSmoothing(TwoDimGrid* grid, double omega, double targetError);
#endif