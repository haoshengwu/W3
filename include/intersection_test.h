#ifndef INTERSECTION_TEST_H
#define INTERSECTION_TEST_H

#include "threedimgridgen.h"
#include "datastructure.h"

/*
!!!PlateCells is coresponding to the final state 3D grid, e.g. grid3d_all_SOL_order, grid3d_all_PFR_order
*/
typedef struct PlateCells
{
  int idx_zone;
  int npcell;
  int nrcell;
  int ntcell;
  //The order is Toroidal-Poloidal-Radial 
  // to ensure the R-P-T sequence
  int*** is_plates;
} PlateCells;


PlateCells* generate_platecells(ThreeDimGrid* grid, int idx_zone);

void set_is_plate_cell(PlateCells* platecells, int ipcell, int ircell, int itcell);

int get_plate_cell(PlateCells* platecells, int ipcell, int ircell, int itcell);

//We assume the grid is the final state, start from outer to inner.
void update_platecells_targets(PlateCells* platecells, ThreeDimGrid* grid);



////We assume the grid is the final state, start from outer to inner.
//The limiter structure is represented by a DDList, it can be CW or CCW (We assume is not closed).
//We also assume the start_phi and end_phi is match with one phi-slice
void update_platecells_limiter(PlateCells* platecells, ThreeDimGrid* grid, DLListNode* limiter, 
                               const double start_phi, const double end_phi);

// n_neu represents the number of neutral cells up to the last cell.
/*
    ----------------------------------------------
    |           Neutral region                   |  
    |                n_rad_neu                   |
    +--------------------------------------------+
    |                                            |
    |           Plasmae region                   |
    |                                            |
    ----------------------------------------------   ^ Radial driection
  Inner                                          Outer
*/
void update_platecells_radneu(PlateCells* platecells, int n_rad_neu);

void free_platecells(PlateCells* platecells);

void write_platecells(PlateCells* platecells, char* filename);



#endif