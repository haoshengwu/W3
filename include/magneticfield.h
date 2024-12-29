#ifndef MAGNETICFIELD_H
#define MAGNETICFIELD_H


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "equilibrium.h"
#include "datastructure.h"
#include "mathbase.h"
//mininum R range used for calculate magnetic field to avoid zero value;
#define MIN_R 0.1 // unit(meter)

typedef struct{
  int nr;
  int nphi;  //currently not use
  int nz;
  double b0r0;
  double dphi;
  double *r;
  double *phi; //currently not use
  double *z;
  //B_rz[i][j][0] is Br, B_rz[i][j][1] is Bz
  double ***Brz;
  //dBrzdx[i][j][0] is dBrdx, dBrzdx[i][j][1] is dBzdx
  double ***dBrzdx;
  //dBrzdy[i][j][0] is dBrdy, dBrzdx[i][j][1] is dBzdy
  double ***dBrzdy;
  //d2Brzdxdy[i][j][0] is d2Brdxdy, d2Brzdxdy[i][j][1] is d2Bzdxdy
  double ***d2Brzdxdy;
} MagFieldTorSys;

/*
follow define a algorithm to select different differential method
*/
typedef void (*diff_2d_fun)(int nx, double *x,  int ny, double *y, double **f, double ***df);

typedef struct {
  const char *name;
  diff_2d_fun func;
} DiffMethodEntry;

diff_2d_fun get_diff_method(const char *name);

/*
follow define a algorithm to select different interpolate method
*/

// 'void *intpl_data' is used to point any intermediate pre-calculated data to speed up the calculation
typedef void (*intpl_2d_fun)(double target_x, double target_y, int nx, double *x,  int ny, double *y, 
                             double ***f, double *value1, double*value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);

typedef struct {
  const char *name;
  intpl_2d_fun func;
} InterpolateMethodEntry;

intpl_2d_fun get_intpl_2d_method(const char *name);

void init_mag_field_torsys(MagFieldTorSys *mag_field);
void free_mag_field_torsys(MagFieldTorSys *mag_field);
void calc_mag_field_torsys(Equilibrium *equ, MagFieldTorSys *mag_field, const char *method);

void get_bt_torsys(MagFieldTorSys *mag_field, const double r0, double *bt);
void get_brz_torsys(MagFieldTorSys*mag_field, const double r, const double z, const char *method, double *br, double *bz);

void write_mag_field_torsys(MagFieldTorSys *mag_field);
void write_brz_torsys(MagFieldTorSys *mag_field);
void write_bphi_torsys(MagFieldTorSys *mag_field);

#endif