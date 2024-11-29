#ifndef ODE_H
#define ODE_H

#include <stdio.h>
#include <math.h>

// define a structure to abstract the funtion f which is the f of dy/dx = f
typedef struct {
  int ndim; // dimension of y, y is y[n].
  void *data; // provide necessary data, not mandatory.
  // a pointer of function used for the calculate of f.
  void (*compute_f)(const int ndim, const double *x, const double *y, double *dydx, void *data);
} ode_function;


// define a structure to abstract the solver which solve the dy/dx = f. 
typedef struct 
{
  double step_size; // step size
  void *solver_data; // Intermediate data
  // calculate next step value y_next
  void (*next_step)(double step_size, const double *x, const double *y, double *y_next, void *solver_data, ode_function* ode_function);
  void (*initialize)(void *solver_data);
  void (*finalize)(void *solver_data);
} ode_solver;


// Following are Specific functions for compute_f

// calcaulte the magnetic B in the R Phi Z coordinates with toroidal sysmmetry. Bpol is calculated by bilinear.
void compute_B_bilinear(const int ndim, const double *x, const double *y, double *dydx, void *data);

// calcaulte the magnetic B in the R Phi Z coordinates with toroidal sysmmetry. Bpol is calculated by 2d cubic Hermite.
void compute_B_hermite2(const int ndim, const double *x, const double *y, double *dydx, void *data);





#endif