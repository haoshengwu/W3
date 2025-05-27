#ifndef TRACERTEST_H
#define TRACERTEST_H
#include "equilibrium.h"




void tracer_test();
void df1dx(const int ndim, const double *x, const double *y, double *dydx, void *data);
double f1(const double x);

void interpolator_test();
double cubic_polynomial_2d(double x, double y);
double quartic_polynomial_2d(double x, double y);
void write_results_2d(int nx, double *x, int ny, double *y, double ***f, const char* filename1, const char* finename2);
void write_results_1d(int nx, double *x, int ny, double *y, double **f, const char* filename);
void line_tracer_test();

// find the estimated value of R with the given psi and Z
double estimate_R_from_z_psi(double r, double z, double psi_value, Equilibrium *equ);
#endif