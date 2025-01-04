#ifndef XPOINT_H
#define XPOINT_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "datastructure.h"
#include "mathbase.h"
#include "equilibrium.h"
#include "calc.h"
#include "magneticfield.h"
#include "input.h"


// Some structs are from DG, the objects with common out are not used.

//Store the info of the extremums (maxmium and minium) for a Xpoint.
//x and y are the point number, t indicates is maximum or minium.
//lvl is the psi value;
struct XPointExtremum {
  int x,y,t;
  double lvl;
};

typedef struct {
  int type,locks;
  //App app;
  int cx1,cy1,cx2,cy2;
  double lvlMin,lvlMax,level,centerX,centerY;
  struct XPointExtremum minMax[4];
  //Group segs;
  //Group gradients[4];
  //int id;
} _XPointInfo;

typedef _XPointInfo* XPointInfo;

//interpl_psi_f can be bound to bilenar_1d/cubicherm_1d/etc. which used to interpolate psi value
typedef void (*interpl_1D_f)(double target_x, double target_y, int nr, double *r,  int nz, double *z,
                      double **psi, double *value, double **dfdx, double **dfdy, double **d2fdxdy);

//interpl_Brz_f can be bound to bilenar_2d/cubicherm_2d/etc. which used to interpolate Br and Bz value.
typedef void (*interpl_2D_f)(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);

void find_xpoint(Equilibrium *equ, int xpoint_number, double **est_xpoint_pos,
                 interpl_1D_f interpl_1D_f, interpl_2D_f interpl_2D_f, _XPointInfo *xpoint);


void test_find_xpoint();
#endif