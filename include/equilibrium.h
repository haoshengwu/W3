#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

//  Define the structure of equlibrium. Currently, it is assuemd the format is geqsk file.
//  
//  Related functions are also clear here.




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
  double rcenter; //R in meter of vacuum toroidal magnetic field BCENTR
  double bcenter; //Vacuum toroidal magnetic field in Tesla at RCENTR
  //O-point
  double rmaxis; //R of magnetic axis in meter
  double zmaxis; //Z of magnetic axis in meter
  double *r, *z;    //  R and Z coordinates, dynamic arrays [m].
  double **psi;      //  Psi value, dynamic arrays [m].

  int Xpoint_num;    //number of X-point, [need to be update]
  double Xpoint_pos[2];

  double minVal;
  double maxVal;

}  Equilibrium;

void init_equilibrium(Equilibrium* equilib); //  initialize 

// read from geqsk file     
void read_equilib_geqdsk(Equilibrium* equilib, const char* geqdsk_file);

// write the equilibrium to dg format [need to be update]
//void write_equilib_dgequ(Equilibrium* equilib, const char* filename);

// find the accurate Xpoint position, est_pos in double array[2] the estimation positions, n i
XPointTest find_Xpoint_from_DG(Equilibrium* equilib, const double* est_pos);

void print_equilibrium(const Equilibrium* equilib);  // print the size of equilibrium

void free_equilibrium(Equilibrium* equilib);         // free dynamic memory

double EqCorrCell(const Equilibrium *equilib,int cx,int cy,double level);

double get_psi_from_rz(const Equilibrium *equilib, double x, double y);

//void find_X_point(const Equilibrium* equilib, const double estimate[2], double accurate[2]); 

//This is to correct the magnetic direction to make sure the magnect field direction is consistent 
//to the current code. Currently, we only support 1 type direction. see magnetic field direction.
//The equlibriums have different coordinate systems that have different positive or negative psi and bcenter.
//This function is correct them. 
//In nutshell, we only support one direction definition, if the equilibirum is not consistent with 
//this defination, we transfer it forcely.
//Current only support lower divertor. Upper divertor need be checked.
void correct_direction_lower_divertor(Equilibrium *equilib);
#endif
