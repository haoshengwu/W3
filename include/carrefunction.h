#ifndef CARREFUNCTION_H
#define CARREFUNCTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datastructure.h"
#include "mathbase.h"
#include "magneticfield.h"

typedef struct
{
/*********************************
* Input parameters
**********************************/
    //point number of the previous curve.
    size_t n_prev_curve;
    //previous curve points prev_curve[n_prev_curve][0] x coordinate, prev_curve[n_prev_curve][0] x coordinate
    double **prev_curve;
    //point number of the curve.
    size_t n_curve;
    //curve points curve[n_curve][0] x coordinate, curve[n_curve][0] x coordinate
    double **curve;
    
    //the number of points for the mesh in X direction
    size_t n_point;
    //The length distributions of the mesh pionts in the previous curve. 
    double *length_prev_points;
    //the points' coordiantes for the mesh which are along the prev_curve.
    double **prev_point_coord;

    //guardlength at the start for curve
    double guard_top;
    //guardlength at the end for curve
    double guard_end;
    //minimum distance between mesh points in the curve
    double pasmin;
/*********************************
* Output parameters
**********************************/
    //The length of the mesh pionts in the curve, is the length from the start point to the point
    double *length_points;
    double **point_coord;
} CarreMeshTube;

void initial_CarreMeshTube(CarreMeshTube *tube, 
                           const size_t n_prev_curve, double **prev_curve, 
                           const size_t n_curve, double **curve, 
                           const size_t n_coordinate, double *distrib_prev_point, double **prev_point_coord,
                           const double guard_top, const double guard_end,const double pasmin, 
                           double *distrib_points, double **point_coordinate);


//CarreOrthoValue is used to store the intermediate values to evaluate the orthogonal property.
typedef struct
{
  size_t n;
  double *ort;
  double *ortpur;
  double *propo;
  double *varr;
  double *tot;
} CarreOrthoValue;
void allocate_CarreOrthoValue(int n, CarreOrthoValue *ortho_value);
void free_CarreOrthoValue(CarreOrthoValue *ortho_value);

/*******************************************************************************
* This function calculate the orthogonalirty of the mesh points on two mesh curves.
* This fuction is refered to Fortran code CARRE/clort.F .
********************************************************************************/
void calc_ortho_CARRE(size_t n_point, double *length_prev_points, double **prev_point_coord,
                      double guard_top, double guard_end, double pasmin,
                      double *length_points, double **point_coord,
                      CarreOrthoValue *ortho_value);

/*******************************************************************************
* This function calculate the mesh points in the curve which have a good orthogonalirty
* using a secant method according to R. Marchand Computer Physics Communications, 1996, 96.2-3: 232-246.
* In carre, the guard length is depend on separatrix. In this code, it dependt on previous cureve.
* Then the guard length are vertorized(different curves have different value). 
********************************************************************************/
void calc_points_CARRE(CarreMeshTube *tube);



/*******************************************************************************
* Other functions from CARRE
********************************************************************************/
int indcrb_CARRE(double **curve, size_t n_curve, double point[2], double d);
double long_CARRE(double **curev, size_t n);
double ruban_CARRE(double **curve, size_t n, double point[2], double d);
void coord_CARRE(double **curve, size_t n, double d, double point[2]);

#endif