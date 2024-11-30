#ifndef ODE_H
#define ODE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datastructure.h"
#include "mathbase.h"
#include "magneticfield.h"

// define a structure to abstract the funtion f which is the f of dy/dx = f
typedef struct {
  //dimension of y, y is y[n].
  int ndim;
  //provide necessary data, not mandatory.
  void *data;
  //rescale the direction and value of dy/dx.
  double *rescale;
  // a pointer of function used for the calculate of f.
  void (*compute_f)(const int ndim, const double *x, const double *y, double *dydx, void *data);
} ode_function;


// define a structure to abstract the solver which solve the dy/dx = f. 
typedef struct 
{
   //step size
  double step_size;
   // Intermediate data
  void *solver_data;
  // calculate next step value y_next
  void (*next_step)(double step_size, const double *x, const double *y, double *y_next, void *solver_data, ode_function *ode_f);
  void (*initialize)(void *solver_data);
  void (*finalize)(void *solver_data);
} ode_solver;


// Following are Specific functions for compute_f

// calcaulte the magnetic B in the R Phi Z coordinates with toroidal sysmmetry. Bpol is calculated by bilinear.
void ode_f_brz_torsys_bilinear(const int ndim, const double *x, const double *y, double *dydx, void *data);

// calcaulte the magnetic B in the R Phi Z coordinates with toroidal sysmmetry. Bpol is calculated by 2d cubic Hermite.
void ode_f_brz_torsys_cubherm(const int ndim, const double *x, const double *y, double *dydx, void *data);

//following are necessary funcstion for euler method
void euler_next_step(double step_size, const double *x, const double *y, double *y_next, void *solver_data, ode_function* ode_function);
void euler_initialize(void *solver_data);
void euler_finalize(void *solver_data);

//following are necessary funcstion for rk45 method
#endif