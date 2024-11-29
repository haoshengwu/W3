#ifndef MATHBASE_H
#define MATHBASE_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define M_PI 3.14159265358979323846


//convertion angle to radians
double deg2rad(double phi);


// df[ix][iy][0] is df/dx, is df[ix][iy][1] is df/dy
void central_2nd_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
void central_4th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
//ToDo
//void central_6th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
//void central_8th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);

void bilenar_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y, 
                double ***f, double *value1, double *value2);
void bicubic_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2);






#endif