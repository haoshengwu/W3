#ifndef THREEDIMGIRDGEN_H
#define THREEDIMGIRDGEN_H

#include "gridzoneinfo.h"
#include "ode.h"
#include "equilibrium.h"
#include "sepdistribution.h"
#include "datastructure.h"
#include "twodimgridgen.h"


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


//THE GridPoint3D IS BASED ON R-Z-Phi coordinate system
typedef struct {
    double r, z, phi;  // Physical coordinates of a Grid point
} GridPoint3D;


// Storage order optimization direction for 3D grid
typedef enum {
    GRID_3D_OPTIMIZE_PRT,  // pol-rad-tor order (ip fastest, then ir, then it)
    GRID_3D_OPTIMIZE_RPT   // rad-pol-tor order (ir fastest, then ip, then it)
} Grid3DOptimization;

// ThreeDimGrid with flexible storage order
typedef struct {
    int npol;       // Logical number of poloidal points
    int nrad;       // Logical number of radial points  
    int ntor;       // Logical number of toroidal points
    
    int cap_npol;   // Allocated number of poloidal points
    int cap_nrad;   // Allocated number of radial points
    int cap_ntor;   // Allocated number of toroidal points
    
    int offset_pol; // Poloidal offset from physical memory origin
    int offset_rad; // Radial offset from physical memory origin
    int offset_tor; // Toroidal offset from physical memory origin
    
    // Grid points stored in 1D array with layout determined by opt_direction
    // For GRID_3D_OPTIMIZE_PRT: points[(offset_tor + it) * cap_npol * cap_nrad + (offset_rad + ir) * cap_npol + (offset_pol + ip)]
    // For GRID_3D_OPTIMIZE_RPT: points[(offset_tor + it) * cap_nrad * cap_npol + (offset_pol + ip) * cap_nrad + (offset_rad + ir)]
    GridPoint3D* points;
    
    Grid3DOptimization opt_direction; // Storage order optimization
    
} ThreeDimGrid;



/*===========================================================================
 * SPECIALIZED CREATION FUNCTIONS
 *===========================================================================*/


 //Storage order: pol-rad-tor (ip fastest, then ir, then it)
 //Suitable for: fixed ir and it, loop over ip

ThreeDimGrid* create_3Dgrid_poloidal_major(int npol, int nrad, int ntor);

//Suitable for: fixed ip and it, loop over ir
//Example use case: radial profile analysis, transport calculations
ThreeDimGrid* create_3Dgrid_radial_major(int npol, int nrad, int ntor);

/**
 * @brief General creation function with optimization direction
 * 
 * This version allocates the grid points array with 32-byte alignment
 * for better cache performance and SIMD operations.
 * 
 * @param npol Number of poloidal points
 * @param nrad Number of radial points
 * @param ntor Number of toroidal points
 * @param opt Storage order optimization (GRID_3D_OPTIMIZE_PRT or GRID_3D_OPTIMIZE_RPT)
 * @return Pointer to new grid or NULL on failure
 */

ThreeDimGrid* create_3Dgrid_optimized_for(int npol, int nrad, int ntor, Grid3DOptimization opt);

/**
 * @brief Free 3D grid and all associated memory
 * 
 * Safely handles NULL pointers and aligned memory allocations.
 * 
 * @param grid Grid to free (can be NULL)
 */
void free_3Dgrid(ThreeDimGrid* grid);

//  ip,ir,it are index
//  Get pointer to grid point
GridPoint3D* get_point_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it);

//  Get r coordinate
double get_r_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it);

//  Get z coordinate
double get_z_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it);

//  Get phi coordinate
double get_phi_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it);

//  Set r z and phi coordinate
void set_point_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double r, double z, double phi);

void set_r_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double r);

void set_z_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double z);

void set_phi_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double phi);

// assign the 2dgrid to the slice of a 3dgrid. it is the toroidal index.
void assign_2D_to_3D_tor_slice(const TwoDimGrid* grid2d, ThreeDimGrid* grid3d, int it, double phim);

//generate the EMC3 (R-Z-Phi CSYS) 3D GRID based on the 2D GRID through magnetic field line tracing.
//The grid2d is one the phim PHI-plane.
//The grid3d is from phi[0] to phi[nphi-1].
void generate_EMC3_3Dgrid_from_2Dgrid_tracing(const TwoDimGrid* grid2d, ThreeDimGrid* grid3d, 
                                              double phim, int nphi, double* phi,  
                                              ode_function* func,ode_solver* solver);


//Write the EMC3 3D grid to XYZ Coordinate System
void write_EMC3_3Dgrid_to_XYZ_CSYS(ThreeDimGrid* grid3d, char* filename);

//Write the EMC3 3D grid to EMC3 required format which is Rad-Pol-Tor order.
//EMC3 3D grid is GRID_3D_OPTIMIZE_RPT rad-pol-tor order (ir fastest, then ip, then it)
void write_EMC3_3Dgrid_to_EMC3_format(ThreeDimGrid* g, char* filename);


//load the EMC3 3D grid (EMC3 format) from the file
//EMC3 3D grid is GRID_3D_OPTIMIZE_RPT rad-pol-tor order (ir fastest, then ip, then it)
ThreeDimGrid* load_EMC3_format_3Dgrid_from_file(char* filename);

#endif















