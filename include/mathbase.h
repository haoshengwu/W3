#ifndef MATHBASE_H
#define MATHBASE_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define M_PI 3.14159265358979323846
#ifndef min
#define min(x,y) ((x)<(y) ? (x) : (y))
#endif
#ifndef max
#define max(x,y) ((x)>(y) ? (x) : (y))
#endif

//convertion angle to radians
double deg2rad(double phi);

/******************************************************
 *    Difference Algorithm
 ******************************************************/
// df[ix][iy][0] is df/dx, is df[ix][iy][1] is df/dy
void central_2nd_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
void central_4th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
//ToDo
//void central_6th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
//void central_8th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);


/******************************************************
 *  Interpolation Algorithm
 *  'void *intpl_data' is used to point any intermediate structure data to speed up the calculation
 *******************************************************/
typedef struct 
{
  int nx;
  int ny;
  double **dfdx;
  double **dfdy;
  double **d2fdxdy;
} bicubic_2d_data;

typedef struct 
{
  double **dfdx;
  double **dfdy;
  double **d2fdxdy;
} CubicHerm2dData;

// 'void *intpl_data' is pointer can be used point to the pre-calculated data to speed up the calculation
void bilenar_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y, 
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);
void bicubic_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);
void cubicherm_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);
#endif