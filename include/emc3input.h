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
void write_tor_axi_sym(ThreeDimGrid* grid, DLListNode* surface, bool reverse, char* filename);

// This function create a non-toroidal Axisymmetric surface file, used to neutral in ADD
// we assume the toroidal range are within the phi range of 3D grid.
void write_non_tor_axi_sym(ThreeDimGrid* grid, DLListNode* surface, bool reverse, 
                           double start_phi, double end_phi, char* filename);

// This function create a toroidal Axisymmetric surface file, used to neutral in ADD
// The range of phi including start, end and delta is decided by the user.
void write_non_tor_axi_sym_usr(DLListNode* surface, bool reverse,
                               double start_phi, double end_phi, double delta_phi, 
                               char* filename);

//we assume surf_start and surf_end have same order: CW or CCW. We don't check this!
void write_one_toroidal_surface(DLListNode* surf_start, double phi_start,
                                DLListNode* surf_end, double phi_end,
                                bool reverse, char* filename);


// Wirte the ip-th (0-base) cell center along radial direction at ir torodial slice to filename
void write_RZ_along_radial_test(ThreeDimGrid* grid, int ip_cell, int it, char* filename);


#endif