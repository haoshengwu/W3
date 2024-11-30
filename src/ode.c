#include "ode.h"

// we assum ndim = 3, y[0] is r (m), y[1] is z(m), y[2] is phi(degree). x is phi. 
void ode_f_brz_torsys_bilinear(const int ndim, const double *x, const double *y, double *dydx, void *data)
{
  MagFieldTorSys *mag_field = (MagFieldTorSys *) data;
  
  double br_tmp, bz_tmp, bt_tmp;
  if( y[0] < MIN_R)
  {
    printf("WARNING: points in the MIN_R: %lf region!\n", MIN_R);
    bt_tmp = mag_field->b0r0 / y[0];
  }
  else
  {
    bt_tmp = mag_field->b0r0 / y[0];
  }
  bilenar_2d(y[0], y[1], mag_field->nr, mag_field->r, mag_field->nz, mag_field->z, 
             mag_field->Brz, &br_tmp, &bz_tmp);
  double dg2rad = M_PI/180.0;
  dydx[0] = y[0] * br_tmp/bt_tmp * dg2rad;
  dydx[1] = y[0] * bz_tmp/bt_tmp * dg2rad;
  dydx[2] = 1.0;
}

void ode_f_brz_torsys_cubherm(const int ndim, const double *x, const double *y, double *dydx, void *data)
{
  //Todo
  return;
}


void euler_next_step(double step_size, const double *x, const double *y, double *y_next, 
                void *solver_data, ode_function* ode_f)
{
  int ndim = ode_f->ndim;
  double dydx[ndim];
  ode_f->compute_f(ndim, x, y, dydx, ode_f->data);
  for (int i=0;i<ndim;i++)
  {
    y_next[i] = y[i] + ode_f->rescale[i] * dydx[i] * step_size;
  }
}

void euler_initialize(void *solver_data)
{
  printf("Euler solver initialized.\n");
}

void euler_finalize(void *solver_data)
{
  printf("Euler solver finalized.\n");
}
