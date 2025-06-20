#ifndef DIVGEO_H
#define DIVGEO_H

#include "gridzone.h"
#include "separatrix.h"
#include "equilibrium.h"
#include "gradpsi.h"
#include "curve.h"
#include "target.h"

typedef struct 
{
  char* name;
  int n_level;
  double* level;
}DGRadRegion;

DGRadRegion* create_DGRadRegion(int n_level, const char* name);
void free_DGRadRegion(DGRadRegion* radregion);

typedef struct 
{
  char* name;
  int n_points;
  double* norm_dist;
}DGPolZone;
DGPolZone* create_DGPolZone(int n_points, const char* name);
void free_DGPolZone(DGPolZone* polzone);


//*.dgo file created by DivGeo
typedef struct
{
  int* n_xpt;
  double** xpt;
  double** opt;
} DivGeoDgo;

//*.str file created by DivGeo\
//NOT USED YET
// typedef struct
// {
//   int n_str;
//   DGClosedStruc** structures;
// } DivGeoStr;

//*.trg file created by DivGeo
//all necessary data for CARRE mesh generation
typedef struct
{
  char* topo;
  int n_target;//total target number
  int* n_target_curve; //array[i] is the total number of points for target_curves[i]
  OldCurve** target_curves;

  // for the radial distribution of magecit surface
  // !!!Unit is Wb/rad
  int n_region;
  DGRadRegion** regions;

  // for the poloidal points distributions
  int n_zone;
  DGPolZone** zones;


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

//n indicate which target curve in trg
TargetDLListCurve* create_target_curve_from_dgtrg(DivGeoTrg* trg, int n);

//From here it is already related with topology
//TOPOLOGY NEED TO BE CONSIDERED
void write_gridzoneinfo_from_dgtrg(DivGeoTrg* trg, Equilibrium* equ, SeparatrixStr* sep);
void write_sn_gridzoneinfo_from_dgtrg(DivGeoTrg* trg, Equilibrium* equ, SeparatrixStr* sep, GradPsiLineStr* gradspsilines);

void free_dgtrg(DivGeoTrg* trg);

//update an exist GridZoneInfo by DivGeoTrg, index means the regions/zons in DivGeoTrg
void update_GridZoneInfo_from_dgtrg(GridZoneInfo* gridzoneinfo, DivGeoTrg* trg, int index);

//update the start points in an exist GridZoneInfo by DivGeoTrg with r and z
void update_GridZoneInfo_start_points(GridZoneInfo* gridzoneinfo, double* r, double* z, int n_point);

void update_GridZoneInfo_end_curve(GridZoneInfo* gridzoneinfo, const TargetDLListCurve* tgt_cur);

//For CORE the start points is also the end point becase of closed surface.
void update_COREGridZoneInfo_end_curve(GridZoneInfo* gridzoneinfo);


void write_polsegms_from_dgtrg(DivGeoTrg* trg, const char* filename);


//update the poloidal normolized points distribution, also the number of point npoint.
// void update_GridZone_pol_norm_distrb(GridZoneInfo* gridzoneinfo, const double* norm_distb, int n_point);

// void update_GridZone_first_pol_points(GridZoneInfo* gridzoneinfo, DLListNode* head);




#endif