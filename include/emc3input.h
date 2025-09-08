#ifndef EMC3INPUT_H
#define EMC3INPUT_H
#include "threedimgridgen.h"
#include "magneticfield.h"
#include "datastructure.h"


//default method to write the file about the BFIELD of a 3d grid.
//cubicherm2d2f method for interpolation
//in EMC3 the order is R-P-T
void write_BFIELD_file_default(ThreeDimGrid* grid, MagFieldTorSys* magnetic, char* name);


//TEST function to write PLATE_MAG file for EMC3,n_neu the number of neutral cells, not the index number.
//THIS FUNCTION IS VERY DANGEROUS, ONLY FRO TEST USE. 
//IN THE FUTURE, A NEW STABLE ONE WILL BE USED.
//grid is the 3D grid, n_neu is the number of cells which are coresponding to neutral
//pol_dir indicate the poloidal direction, 
//it can be from inner to outer (W3 order) or from outer to inner (TLunt order),
//THIS ORDER MUST BE CONSISTENT WITH THE GRID, HOWEVER, WE DON'T CHECK THAT.
//SO THIS FUNCTION IS DANGEROUS. 
//In the future, a new statble algorithm will be developed.
/*

    ----------------------------------------------
    |           Neutral region                   |  
    |                n_neu                       |
    +--------------------------------------------+
    |                                            |
    |           Plasmae region                   |
    |                                            |
    ----------------------------------------------   ^ Radial driection
  Inner                                          Outer

    pol_dir==0 is from Inner to Outer ===>
    pol_dir==1 is from Outer to Inner <===
*/
void write_PLATE_MAG_file_test(ThreeDimGrid* grid, int n_neutral_cells, int pol_dir, int idx_zone, char* filename);



// This function create a toroidal Axisymmetric surface file, used to neutral in ADD
void write_axis_sys_surface_default(ThreeDimGrid* grid, DLListNode* surface, bool reverse, char* filename);


#endif