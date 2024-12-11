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
             mag_field->Brz, &br_tmp, &bz_tmp, NULL, NULL, NULL);

  //printf("debug br_tmp: %lf, bz_tmp: %lf\n", br_tmp, bt_tmp);
  double dg2rad = M_PI/180.0;
  dydx[0] = y[0] * br_tmp/bt_tmp * dg2rad;
  dydx[1] = y[0] * bz_tmp/bt_tmp * dg2rad;
  if(ndim ==3)
  {
    dydx[2] = 1.0;
  }
}

void ode_f_brz_torsys_bicubic(const int ndim, const double *x, const double *y, double *dydx, void *data)
{
  //Todo
  return;
}

void ode_f_brz_torsys_cubicherm(const int ndim, const double *x, const double *y, double *dydx, void *data)
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
  cubicherm_2d(y[0], y[1], mag_field->nr, mag_field->r, mag_field->nz, mag_field->z, 
             mag_field->Brz, &br_tmp, &bz_tmp, mag_field->dBrzdx, mag_field->dBrzdy, mag_field->d2Brzdxdy);
  double dg2rad = M_PI/180.0;
  //printf("debug br_tmp: %lf, bz_tmp: %lf\n", br_tmp, bt_tmp);
  dydx[0] = y[0] * br_tmp/bt_tmp * dg2rad;
  dydx[1] = y[0] * bz_tmp/bt_tmp * dg2rad;
  if(ndim ==3)
  {
    dydx[2] = 1.0;
  }

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

//Butcher’s (1964) fifth-order RK method
void brk5_initialize(void *solver_data)
{
  RKSolverData *brk5_data = (RKSolverData *) solver_data;
  brk5_data->order = 5;
  brk5_data->stages = 6;
  brk5_data->A = allocate_2d_array(6,6);
  for(int i=0; i<6; i++)
  {
    for(int j=0; j<6; j++)
    {
      brk5_data->A[i][j] = 0.0;
    }
  }
  brk5_data->B = calloc(6, sizeof(double));
  brk5_data->B_ALT = NULL;
  brk5_data->C = calloc(6, sizeof(double));
  brk5_data->current_step_size = 0.0;
  brk5_data->tollor=0.0;
  brk5_data->A[1][0] = 1.0 / 4.0;
  brk5_data->A[2][0] = 1.0 / 8.0; brk5_data->A[2][1] = 1.0 / 8.0;
  brk5_data->A[3][1] = -1.0/2.0; brk5_data->A[3][2] = 1.0;
  brk5_data->A[4][0] = 3.0/16.0; brk5_data->A[4][3] = 9.0/16.0;
  brk5_data->A[5][0] = -3.0/7.0; brk5_data->A[5][1] = 2.0/7.0;
  brk5_data->A[5][2] = 12.0/7.0; brk5_data->A[5][3] =-12.0/7.0;
  brk5_data->A[5][4] = 8.0/7.0;   
  
  brk5_data->C[0] = 0.0;
  brk5_data->C[1] = 1.0 / 4.0;
  brk5_data->C[2] = 1.0 / 4.0;
  brk5_data->C[3] = 1.0 / 2.0;
  brk5_data->C[4] = 3.0 / 4.0;
  brk5_data->C[5] = 1.0;

  brk5_data->B[0] = 7.0 / 90.0;
  brk5_data->B[1] = 0.0;
  brk5_data->B[2] = 32.0 / 90.0;
  brk5_data->B[3] = 12.0 / 90.0;
  brk5_data->B[4] = 32.0 / 90.0;
  brk5_data->B[5] = 7.0 / 90.0;
  printf("BRK5 Solver data Initialize.\n");
}

void brk5_finalize(void *solver_data)
{
  RKSolverData *brk5_data = (RKSolverData *) solver_data;
  free(brk5_data->B);
  free(brk5_data->C);
  free_2d_array(brk5_data->A);
  printf("BRK5 Solver data finalized.\n");
}

void brk5_next_step(double step_size, const double *x, const double *y, double *y_next,
                    void *solver_data, ode_function* ode_f)
{
  RKSolverData *brk5_data = (RKSolverData *) solver_data;
  int ndim = ode_f->ndim;
  int order = brk5_data->order;
  int stages = brk5_data->stages;

  // the dydx for 6 stages of each dimenstions.
  double dydx[6][ndim];
  // temperary y for 6 stages of each dimenstions.
  double y_tmp[ndim];
  double x_tmp;
 
  //step 1, k1
  ode_f->compute_f(ndim, x, y, dydx[0], ode_f->data);
  for (int i=0; i<ndim; i++)
  {
    dydx[0][i] = dydx[0][i] * ode_f->rescale[i];
  }
  //printf("debug: dydx[0] calculated.\n");


  //step 2, k2
  x_tmp = *x + brk5_data->C[1] * step_size;
  for (int i = 0; i < ndim; i++) 
  {
      y_tmp[i] = y[i] + step_size * brk5_data->A[1][0] * dydx[0][i];
  }
  ode_f->compute_f(ndim, &x_tmp, y_tmp, dydx[1], ode_f->data);
  for (int i=0; i<ndim; i++)
  {
    dydx[1][i] = dydx[1][i] * ode_f->rescale[i];
  }
  //printf("debug: dydx[1] calculated.\n");

  // Step 3: k3
  x_tmp = *x + brk5_data->C[2] * step_size;
  for (int i = 0; i < ndim; i++) 
  {
    y_tmp[i] = y[i] + step_size * (brk5_data->A[2][0] * dydx[0][i] + brk5_data->A[2][1] * dydx[1][i]);
  }
  ode_f->compute_f(ndim, &x_tmp, y_tmp, dydx[2], ode_f->data);
  for (int i=0; i<ndim; i++)
  {
    dydx[2][i] = dydx[2][i] * ode_f->rescale[i];
  }
  //printf("debug: dydx[3] calculated.\n");

  // Step 4: k4
  x_tmp = *x + brk5_data->C[3] * step_size;
  for (int i = 0; i < ndim; i++) 
  {
    y_tmp[i] = y[i] + step_size * (brk5_data->A[3][1] * dydx[1][i] + brk5_data->A[3][2] * dydx[2][i]);
  }
  ode_f->compute_f(ndim, &x_tmp, y_tmp, dydx[3], ode_f->data);
  for (int i=0; i<ndim; i++)
  {
    dydx[3][i] = dydx[3][i] * ode_f->rescale[i];
  }
  //printf("debug: dydx[4] calculated.\n");

  // Step 5: k5
  x_tmp = *x + brk5_data->C[4] * step_size;
  for (int i = 0; i < ndim; i++) 
  {
    y_tmp[i] = y[i] + step_size * (brk5_data->A[4][0] * dydx[0][i] + brk5_data->A[4][3] * dydx[3][i]);
  }
  ode_f->compute_f(ndim, &x_tmp, y_tmp, dydx[4], ode_f->data);

  for (int i=0; i<ndim; i++)
  {
    dydx[4][i] = dydx[4][i] * ode_f->rescale[i];
  }
  //printf("debug: dydx[5] calculated.\n");

  x_tmp = *x + brk5_data->C[5] * step_size;
  for (int i = 0; i < ndim; i++) 
  {
    y_tmp[i] = y[i] + step_size * (
            brk5_data->A[5][0] * dydx[0][i] +
            brk5_data->A[5][1] * dydx[1][i] +
            brk5_data->A[5][2] * dydx[2][i] +
            brk5_data->A[5][3] * dydx[3][i] +
            brk5_data->A[5][4] * dydx[4][i]);
  }
  ode_f->compute_f(ndim, &x_tmp, y_tmp, dydx[5], ode_f->data);
  
  for (int i=0; i<ndim; i++)
  {
    dydx[5][i] = dydx[5][i] * ode_f->rescale[i];
  }
  //printf("debug: dydx[6] calculated.\n");


  for (int i = 0; i < ndim; i++) 
  {
    y_next[i] = y[i] + step_size * (
        brk5_data->B[0] * dydx[0][i] +
        brk5_data->B[1] * dydx[1][i] +
        brk5_data->B[2] * dydx[2][i] +
        brk5_data->B[3] * dydx[3][i] +
        brk5_data->B[4] * dydx[4][i] +
        brk5_data->B[5] * dydx[5][i]);
  }

}