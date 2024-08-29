#ifndef __EQUILIBRIUM_H_
#define __EQUILIBRIUM_H_

//  Define the structure of equlibrium. Currently, it is assuemd the format is geqsk file.
//  
//  Related functions are also clear here.

#include <stdlib.h>
#include <string.h>


// There are some codes directly from DivGeo, the reason is I will to use benefits from DivGeo
// It have a strong capaticy about object oriented programming (OOP). In the future, this project
// will be transfered to C++.


// Some structs are from DG, the objects with common out are not used.
struct _XPointMinMax {
  int x,y,t;
  double lvl;
};

// From DG
struct _XPointTest {
  int type,locks;
  //App app;
  int cx1,cy1,cx2,cy2;
  double lvlMin,lvlMax,level,centerX,centerY;
  struct _XPointMinMax minMax[4];
  //Group segs;
  //Group gradients[4];
  //int id;
};

// From DG
typedef struct _XPointTest* XPointTest;

typedef struct {

  int  nw;    //  number of points along the R coordinate . 
  int  nh;    //  number of points along the Z coordinate .
  double simag;  //poloidal flux at magnetic axis. [Weber/rad]
  double sibry;  //poloidal flux at the plasma boundary [Weber/rad].
  double *r, *z;    //  R and Z coordinates, dynamic arrays [m].
  double **psi;      //  Psi value, dynamic arrays [m].

  int Xpoint_num;    //number of X-point, [need to be update]
  double Xpoint_pos[2];

  int Opoint_num;    //number of X-point, [need to be update]
  double Opoint_pos[2];

  double minVal;
  double maxVal;

}  Equilibrium;

void init_equilibrium(Equilibrium* equilib); //  initialize 

// read from geqsk file     
void read_equilib_geqdsk(Equilibrium* equilib, const char* geqdsk_file);

// write the equilibrium to dg format [need to be update]
//void write_equilib_dgequ(Equilibrium* equilib, const char* filename);

// find the accurate Xpoint position, est_pos in double array[2] the estimation positions, n i
void find_Xpoint(Equilibrium* equilib, const double* est_pos);

void print_equilibrium(const Equilibrium* equilib);  // print the size of equilibrium

void free_equilibrium(Equilibrium* equilib);         // free dynamic memory

double EqCorrCell(const Equilibrium *equilib,int cx,int cy,double level);
//void find_X_point(const Equilibrium* equilib, const double estimate[2], double accurate[2]); 
#endif
