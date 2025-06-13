#ifndef DIVGEO_H
#define DIVGEO_H

#include "structuredgrid.h"
#include "separatrix.h"
#include "equilibrium.h"
#include "gradpsi.h"
#include "curve.h"
#include "structuredgrid.h"

typedef struct 
{
  int n_point;
  double** points;

} DGClosedStruc;

DGClosedStruc* create_DGClosedStruc();
void free_DGClosedStruc(DGClosedStruc* struc);

typedef struct 
{
  char* name;
  int n_level;
  double* level;
}DGRegion;

DGRegion* create_DGRegion(int n_level, const char* name);
void free_DGRegion(DGRegion* region);

typedef struct 
{
  char* name;
  int n_points;
  double* norm_dist;
}DGZone;
DGZone* create_DGZone(int n_points, const char* name);
void free_DGZone(DGZone* zone);


//*.dgo file created by DivGeo
typedef struct
{
  int* n_xpt;
  double** xpt;
  double** opt;
} DivGeoDgo;

//*.str file created by DivGeo
typedef struct
{
  int n_str;
  DGClosedStruc** structures;
} DivGeoStr;

//*.trg file created by DivGeo
//all necessary data for CARRE mesh generation
typedef struct
{
  char* topo;
  int n_target;//total target number
  int* n_target_curve; //array[i] is the total number of points for target_curves[i]
  Curve** target_curves;

  // for the poloidal points distributions 
  int n_region;
  DGRegion** regions;

  // for the radial distribution of magecit surface
  // !!!Unit is Wb/rad
  int n_zone;
  DGZone** zones;


  double* dltr1;
  double* dltrn;
  int* npr;
  
  double pntrat;

  double* dltp1;
  double* dltpn;
  int* nptseg;
} DivGeoTrg;

DivGeoTrg* create_dgtrg(void);
int load_dgtrg_from_file(DivGeoTrg* trg, const char* filename);

//From here it is already related with topology
//TOPOLOGY NEED TO BE CONSIDERED
void write_dgtrg_to_input(DivGeoTrg* trg, Equilibrium* equ, SeparatrixStr* sep);
void write_dgtrg_to_sn_input(DivGeoTrg* trg, Equilibrium* equ, SeparatrixStr* sep, GradPsiLineStr* gradspsilines);

void free_dgtrg(DivGeoTrg* trg);



//update an exist GridZone by DivGeoTrg, index means the regions/zons in DivGeoTrg
void update_GridZone_from_dgtrg(GridZone* gridzone, DivGeoTrg* trg, int index);

//update the start points in an exist GridZone by DivGeoTrg with r and z
void update_GridZone_start_points(GridZone* gridzone, double* r, double* z, int n_points);




#endif