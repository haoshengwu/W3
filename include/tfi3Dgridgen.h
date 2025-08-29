#ifndef TFT3DGRIDGEN
#define TFT3DGRIDGEN


#include "threedimgridgen.h"

//Generate the 3D grid for neutral transport based on the TFI
//grid3d_neu is the output of neutral 3d grid
//INPUT:
//grid3d_plasma is the 3D plasma grid 
//inner_tgt and outer_tgt are the DDList target,the values are the same with the target_curves in DivGeoTrg
// n_neu_distrb is number for neutral grid in radial direction
// neu_distrb is the distribution that neu_distrb[0]=0.0 and neu_distrb[n_neu_distrb-1]=1.0
// Here we assume that neu_bnd should be connected with the head or tail of inner_tgt and outer_tgt
void generate_EMC3_neutral_3Dgrid_TFI(ThreeDimGrid* grid3d_neu,
                                      const ThreeDimGrid* grid3d_plasma,
                                      DLListNode* inner_tgt, DLListNode* outer_tgt, DLListNode* neu_bnd,
                                      const int n_neu_distrb, const double* neu_distrb,
                                      ode_function* func, ode_solver* solver);

#endif