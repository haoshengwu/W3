#ifndef TRACERTEST_H
#define TRACERTEST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "equilibrium.h"
#include "magneticsurface.h"
#include "datastructure.h"
#include "linetrace.h"
#include "magneticfield.h"
#include "ode.h"
#include "basemesh.h"
#include "mathbase.h"


void tracer_test();
void df1dx(const int ndim, const double *x, const double *y, double *dydx, void *data);
double f1(const double x);

void interpolator_test();
double cubic_polynomial_2d(double x, double y);
double quartic_polynomial_2d(double x, double y);
void write_results_2d(int nx, double *x, int ny, double *y, double ***f, const char* filename1, const char* finename2);
void write_results_1d(int nx, double *x, int ny, double *y, double **f, const char* filename);

#endif