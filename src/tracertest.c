/* This code is just test the magnetic field line tracer.
 * Should not be included in source code. 
*/
#include "tracertest.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "magneticsurface.h"
#include "datastructure.h"
#include "linetrace.h"
#include "magneticfield.h"
#include "ode.h"
#include "carrefunction.h"
#include "mathbase.h"
/****************************************************************************
 * Some static functions used for test the tracer.
****************************************************************************/
void df1dx(const int ndim, const double *x, const double *y, double *dydx, void *data)
{
    if(ndim != 1)
    {
        printf("This is test function only support ndim = 1! \n");
        assert(ndim ==1);
    }
    *dydx = -2*pow(*x,3.0) + 12*pow(*x,2.0) - 20 * (*x) + 8.5;
};

double f1(const double x)
{
    //y = -0.5x^4 + 4x^3 - 10x^2 + 8.5x + 1
    double y = -0.5*pow(x,4.0) + 4*pow(x,3.0) - 10 * pow(x,2.0) + 8.5 *(x) + 1;
    return y;
}


void tracer_test()
{
  printf("*******************************************\n");
  printf("Begint to test the tracer\n");

/*****************************************************
test euler method and brk5 tracing for 
function f1 that y=f1(x)
                 y = -0.5x^4 + 4x^3 - 10x^2 + 8.5x + 1
           thus  dy/dx = -2x^3 + 12x^2 - 20x + 8.5
******************************************************/
  double size = 0.1;
  int step = 40;
  double direction = 1.0;


  ode_function ode_df1dx = {
    .ndim = 1,
    .data = NULL,
    .rescale = &direction,
    .compute_f = df1dx
  };

  ode_solver euler_solver = {
    .step_size = size,
    .solver_data = NULL,
    .next_step = euler_next_step,
    .initialize = euler_initialize,
    .finalize = euler_finalize,
  };

 
  double *x1 = (double *)malloc((step+1) * sizeof(double));
  x1[0]=0.0;
  double **y_euler = allocate_2d_array(step+1,1);
  y_euler[0][0] = 1.0;

  printf("debug in tracer\n");
  for(int i=1;i<step+1;i++)
  {
    //printf("debug in loop, i: %d\n", i);
    x1[i] = x1[i-1] + euler_solver.step_size;
    //printf("debug in loop, x1[i]: %f\n", x1[i-1]);
    //printf("debug in loop, y_euler[i]: %f\n", y_euler[i-1][0]);
    euler_solver.next_step(euler_solver.step_size, &x1[i-1], y_euler[i-1], y_euler[i], NULL, &ode_df1dx);
  };
  printf("debug in tracer\n");
  
  RKSolverData brk45_data;

  ode_solver brk45_solver =
  {
    .step_size = size,
    .solver_data = NULL,
    .next_step = brk5_next_step,
    .initialize = brk5_initialize,
    .finalize = brk5_finalize
  };
  brk45_solver.initialize(&brk45_data);

  double **y_brk5 = allocate_2d_array(step+1,1);
  y_brk5[0][0] = 1.0;

  printf("debug in tracer\n");
  for(int i=1;i<step+1;i++)
  {
    //printf("debug in loop, i: %d\n", i);
    brk45_solver.next_step(brk45_solver.step_size, &x1[i-1], y_brk5[i-1], y_brk5[i], &brk45_data, &ode_df1dx);
  };
  printf("debug in tracer\n");



//out put euler method result
  const char *filename1="euler_tracer";
  FILE* file1 = fopen(filename1, "w");
  fprintf(file1, "#x, y_euler\n");

  for (int i=0; i<step+1; i++)
  {
    fprintf(file1, "%.12f %.12f\n", x1[i],y_euler[i][0]);
  }
  fclose(file1);
  printf("write the tracing line in %s\n", filename1);

//out put brk5 method result
  const char *filename2="brk5_tracer";
  FILE* file2 = fopen(filename2, "w");
  fprintf(file2, "#x, y_brk5\n");

  for (int i=0; i<step+1; i++)
  {

    fprintf(file2, "%.12f %.12f\n", x1[i],y_brk5[i][0]);
  }
  fclose(file2);
  printf("write the tracing line in %s\n", filename2);



//out put accurate result
  const char *filename3="accurate_tracer";
  FILE* file3 = fopen(filename3, "w");
  fprintf(file3, "#x, y\n");

  double delta_x_accurate = 0.01;
  for (int i=0; i<400; i++)
  {
    double y;
    y=f1(i*delta_x_accurate);
    fprintf(file3, "%.12f %.12f\n", i*delta_x_accurate,y);
  }
  fclose(file3);
  printf("write the tracing line in %s\n", filename3);
  
  const char *filename4="tracer_error";
  FILE* file4 = fopen(filename4, "w");
  fprintf(file4, "#x, y, euler, euler_abs, euler_rel, brk5, brk5_abs, brk5_rel\n");
  for (int i=0; i<step+1; i++)
  {
    double y;
    y=f1(x1[i]);
    double euler_abs, euler_rel, brk5_abs, brk5_rel;
    euler_abs = fabs(y_euler[i][0]-y);
    euler_rel = euler_abs/fabs(y);
    brk5_abs =fabs(y_brk5[i][0]-y);
    brk5_rel = brk5_abs/fabs(y);
    fprintf(file4, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", 
            x1[i], y, y_euler[i][0], euler_abs, euler_rel, y_brk5[i][0],brk5_abs, brk5_rel);
  }
  fclose(file4);
  printf("write the error of tracing line in %s\n", filename4);
  
  
  //   free(x2);
//   free_2d_array(line2);
  brk45_solver.finalize(&brk45_data);
  free(x1);
  free_2d_array(y_euler);
  free_2d_array(y_brk5);
  printf("Finish test the tracer\n");
  printf("*******************************************\n");
}



double cubic_polynomial_2d(double x, double y)
{
   return 3.0 * pow(x,3.0) - 2.0 * pow(y,3.0) + 1.5 * pow(x,2.0) * y - 0.5 * x * pow(y,2.0)+ 2.0 * x + 3.0 * y + 5.0;
}

double quartic_polynomial_2d(double x, double y)
{
  return 2.5 * pow(x,4.0) + 1.8 * pow(y,4.0) - 1.2 * pow(x,2.0) * pow(y,2.0) + 0.8 * x * y + 4.0 * x + 2.0 * y + 10.0;
}

void interpolator_test()
{
  
  printf("*******************************************\n");
  printf("Begin test interpolator\n");

  int nx = 11;
  int ny = 21;
  double x[nx];
  double y[ny];
  double dx = 0.01;
  double dy = 0.005;

  double ***f = allocate_3d_array(nx,ny,2);

/*********************************************
 * basic matrix for interpolation            *
**********************************************/

  for(int i=0; i<nx; i++)
  {
    double x0 = 0.0;
    x[i] = dx * i + x0;
  }

  for(int j=0; j<ny; j++)
  {
    double y0 = -0.5; 
    y[j] = dy * j + y0;
  }

  for(int i=0; i< nx; i++)
  {
    for(int j=0; j< ny; j++)
    {
      f[i][j][0] = cubic_polynomial_2d(x[i],y[j]);
      f[i][j][1] = quartic_polynomial_2d(x[i],y[j]);
    }
  }
  write_results_2d(nx, x, ny, y, f, "origin1", "origin2");
  printf("finish write origin resutls\n");
/*********************************************
 * Interpolated matrix                       *
**********************************************/
  int nx1 = 31;
  int ny1 = 61;
  double x1[nx1];
  double y1[ny1];
  double dx1 = (x[nx-2] - x[1]) / (nx1 - 1);
  double dy1 = (y[ny-2] - y[1]) / (ny1 - 1);
  printf("dx1: %f, dy1: %f\n", dx1, dy1);

  for (int i=0; i<nx1; i++)
  {
    x1[i] = x[1] + dx1 * i;
  }

  for (int j=0; j<ny1; j++)
  {
    y1[j] = y[1] + dy1 * j;
    //printf("y1[j]: %f\n", y1[j]);
  }
  x1[0] = x[1];
  y1[0] = y[1];

  x1[nx1-1] = x[nx-2];
  y1[ny1-1] = y[ny-2];

  double ***f_bilinear=allocate_3d_array(nx1,ny1,2);
  double ***f_cubicherm=allocate_3d_array(nx1,ny1,2);

  for(int i=0; i<nx1; i++)
  {
    for(int j=0; j<ny1; j++)
    {
      // printf("i: %d, j: %d\n", i, j);
      // printf("x1[i]: %f, y1[j]: %f\n", x1[i], y1[j]);
      double *f0 = &f_cubicherm[i][j][0];
      double *f1 = &f_cubicherm[i][j][1];
      cubicherm2d2f(x1[i],y1[j],nx, x, ny, y, f, f0, f1, NULL, NULL, NULL);
      // printf("finish cubicherm2d2f interpolator\n");
      f0 = &f_bilinear[i][j][0];
      f1 = &f_bilinear[i][j][1];
      bilenar2d2f(x1[i],y1[j],nx, x, ny, y, f, f0, f1, NULL, NULL, NULL);
      // printf("finish bilenar2d2f interpolator\n");
    }
  }

/*********************************************
* Accurate matrix                            *
**********************************************/
  double ***f_accurate=allocate_3d_array(nx1,ny1,2);
  for(int i=0; i<nx1; i++)
  {
    for(int j=0; j<ny1; j++)
    {
      f_accurate[i][j][0] = cubic_polynomial_2d(x1[i],y1[j]);
      f_accurate[i][j][1] = quartic_polynomial_2d(x1[i],y1[j]);
    }
  }  
  
  write_results_2d(nx1, x1, ny1, y1, f_bilinear, "bilinear1", "bilinear2");
  printf("finish write bilinear resutls\n");

  write_results_2d(nx1, x1, ny1, y1, f_cubicherm, "cubicherm1", "cubicherm2");
  printf("finish write cubicherm resutls\n");

  write_results_2d(nx1, x1, ny1, y1, f_cubicherm, "accurate1", "accurate2");
  printf("finish write accurate resutls\n");

  double ***error_abs=allocate_3d_array(nx1,ny1,2);
  double ***error_rel=allocate_3d_array(nx1,ny1,2);

  for(int i=0; i<nx1; i++)
  {
    for(int j=0; j<ny1; j++)
    {
      error_abs[i][j][0] = fabs(f_bilinear[i][j][0]-f_accurate[i][j][0]);
      error_abs[i][j][1] = fabs(f_bilinear[i][j][1]-f_accurate[i][j][1]);
      error_rel[i][j][0] = fabs(f_bilinear[i][j][0]-f_accurate[i][j][0])/fabs(f_accurate[i][j][0]);
      error_rel[i][j][1] = fabs(f_bilinear[i][j][1]-f_accurate[i][j][1])/fabs(f_accurate[i][j][1]);
    }
  }

  write_results_2d(nx1, x1, ny1, y1, error_abs, "bi_error_abs1", "bi_error_abs2");
  write_results_2d(nx1, x1, ny1, y1, error_rel, "bi_error_rel1", "bi_error_rel2");
  printf("finish write bilinear error calculation\n");



  for(int i=0; i<nx1; i++)
  {
    for(int j=0; j<ny1; j++)
    {
      error_abs[i][j][0] = fabs(f_cubicherm[i][j][0]-f_accurate[i][j][0]);
      error_abs[i][j][1] = fabs(f_cubicherm[i][j][1]-f_accurate[i][j][1]);
      error_rel[i][j][0] = fabs(f_cubicherm[i][j][0]-f_accurate[i][j][0])/fabs(f_accurate[i][j][0]);
      error_rel[i][j][1] = fabs(f_cubicherm[i][j][1]-f_accurate[i][j][1])/fabs(f_accurate[i][j][1]);
    }  
  }
  write_results_2d(nx1, x1, ny1, y1, error_abs, "cu_error_abs1", "cu_error_abs2");
  write_results_2d(nx1, x1, ny1, y1, error_rel, "cu_error_rel1", "cu_error_rel2");
  printf("finish write cubicherm error calculation\n");

  free_3d_array(f);
  free_3d_array(f_bilinear);
  free_3d_array(f_cubicherm);
  free_3d_array(f_accurate);
  free_3d_array(error_abs);
  free_3d_array(error_rel);

  printf("Finish test interpolator\n");
  printf("*******************************************\n");

  return;
}

void write_results_2d(int nx, double *x, int ny, double *y, double ***f, const char* filename1, const char* filename2)
{
  FILE* file1 = fopen(filename1, "w");
  fprintf(file1, "#x, y, value1\n");
  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      fprintf(file1, "%.15f %.15f %.15f\n", x[i], y[j], f[i][j][0]);
    }
  }
  fclose(file1);

  FILE* file2 = fopen(filename2, "w");
  fprintf(file2, "#x, y, value2\n");
  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      fprintf(file2, "%.15f %.15f %.15f\n", x[i], y[j], f[i][j][1]);
    }
  }
  fclose(file2);
}


void write_results_1d(int nx, double *x, int ny, double *y, double **f, const char* filename)
{
  FILE* file = fopen(filename, "w");
  fprintf(file, "#x, y, value1\n");
  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      fprintf(file, "%.15f %.15f %.15f\n", x[i], y[j], f[i][j]);
    }
  }
  fclose(file);
}

void line_tracer_test()
{
  InputPara w3_input;
  init_inputpara(&w3_input);
  print_inputpara(&w3_input);
  Equilibrium dtt_example;
  
  init_equilibrium(&dtt_example);
 
  read_equilib_geqdsk(&dtt_example,w3_input.equilibrium_file);
  
  print_equilibrium(&dtt_example);
  XPointTest xp;
  xp = find_Xpoint(&dtt_example, w3_input.xpt_estimation);

  double value;
  double x_point = 2.75;
  double y_point = 0;
  value = get_psi_from_rz(&dtt_example,x_point, y_point );
  printf("psi at %lf %lf is %lf\n",x_point, y_point, value);

  // x_point = 3.0;
  // y_point = 0;
  // value = get_psi_from_rz(&dtt_example,x_point, y_point );
  // printf("psi at %lf %lf is %lf\n",x_point, y_point, value);

//**********Test four lines of the separatrix tracing********************/

  DLListNode* line_list[4]={NULL};
  for (int i=0; i<4; i++)
  {
    cal_separatrix_line(&dtt_example, xp, i, &(line_list[i]));
  }
  write_DDList(line_list[0],"sep0");
  write_DDList(line_list[1],"sep1");
  write_DDList(line_list[2],"sep2");
  write_DDList(line_list[3],"sep3");

  for (int i=0; i<4; i++)
  {
    free_DLList(line_list[i]);
  }



  printf("*******************************************\n");
  printf("Begin test magnetic field line tracer\n");


/*test magnetic field line calculation*/
  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);
  write_mag_field_torsys(&test_magfield);

/*****************************************************
* Define the magnetic and ODE solver
******************************************************/
  double direction[3]={1.0,1.0,1.0};
  RKSolverData brk45_data;

  double stepsize = 0.1;

  ode_function ode_func = {
    .ndim = 3,
    .data = &test_magfield,
    .rescale = direction,
    .compute_f = ode_f_brz_torsys_cubicherm,
  };

  ode_solver brk45_solver =
  {
    .step_size = stepsize,
    .solver_data = &brk45_data,
    .next_step = brk5_next_step,
    .initialize = brk5_initialize,
    .finalize = brk5_finalize
  };

    ode_solver euler_solver =
  {
    .step_size = stepsize,
    .solver_data = NULL,
    .next_step = euler_next_step,
    .initialize = euler_initialize,
    .finalize = euler_finalize
  };

  brk45_solver.initialize(&brk45_data);


/*****************************************************
* Define necessary input for tracer
******************************************************/


  double start_x = 0.0;
  double prev_x;
  double x;

  prev_x = start_x;
  x = start_x;

  int turn = 1000;
  int one_turn = (int)round(360.0 / stepsize);
  printf("points in one turn: %d\n", one_turn);
  double start_point[3] ={2.6, 0.2, 0.0};

  double *prev_point=calloc(3, sizeof(double));
  double *point=calloc(3, sizeof(double));

  prev_point[0] = start_point[0];
  prev_point[1] = start_point[1];
  prev_point[2] = start_point[2];

  point[0] = start_point[0];
  point[1] = start_point[1];
  point[2] = start_point[2];

/*****************************************************
* Define output of the tracer
******************************************************/
  double **tracer_result = allocate_2d_array(turn,3);


/*****************************************************
* Begin to trace and writ the results
******************************************************/
  const char *filename_debug= "line_tracer_test_debug";
  FILE* file_debug = fopen(filename_debug, "w");
  fprintf(file_debug, "#R,Z\n");
  fprintf(file_debug, "%.15f %.15f\n", point[0], point[1]);

  for (int j = 0; j<turn; j++)
  {
    for (int i=0; i<one_turn; i++)
    {
      x = x + brk45_solver.step_size;
      brk45_solver.next_step(brk45_solver.step_size, &prev_x, prev_point, point, &brk45_data, &ode_func);
      //euler_solver.next_step(euler_solver.step_size, &prev_x, prev_point, point, NULL, &ode_func);
      prev_x = x;
      prev_point[0] = point[0];
      prev_point[1] = point[1];
      prev_point[2] = point[2];
      fprintf(file_debug, "%.15f %.15f\n", point[0], point[1]);
    }
    printf("debug Turn: %d\n", j);

    x = x - 360.0;
    prev_x = prev_x - 360.0;
    tracer_result[j][0] = point[0];
    tracer_result[j][1] = point[1];
    cubicherm2d1f(point[0], point[1], dtt_example.nw, dtt_example.r, dtt_example.nh, dtt_example.z,
                 dtt_example.psi, &(tracer_result[j][2]), NULL, NULL, NULL);
  }
  
  fclose(file_debug);

  double psi_start;
  cubicherm2d1f(start_point[0], start_point[1], dtt_example.nw, dtt_example.r, dtt_example.nh, dtt_example.z,
                 dtt_example.psi, &(psi_start), NULL, NULL, NULL);
  

  const char *filename= "line_tracer_test_results";
  FILE* file = fopen(filename, "w");
  fprintf(file, "#R, Z, Psi\n");
  for (int i=0; i<turn; i++)
  {
    fprintf(file, "%.15f %.15f %.15f\n", tracer_result[i][0], tracer_result[i][1], tracer_result[i][2]);
  }
  fclose(file);
  printf("Finish write the line tracer results\n");


  const char *filename1= "line_tracer_test_error";
  FILE* file1 = fopen(filename1, "w");
  fprintf(file1, "#turn, rel delta psi , delta R\n");
  for (int i=0; i<turn; i++)
  {
    double r;
    r = estimate_R_from_z_psi(tracer_result[i][2], start_point[0], start_point[1], &dtt_example);
    //fprintf(file1, "%.d %.15f \n", i+1, r);
    fprintf(file1, "%.d %.15f %.15f\n", i+1, fabs(tracer_result[i][2]-psi_start)/psi_start, fabs(r-start_point[0])*100);
  }
  fclose(file1);
  printf("Finish write the error of line tracer results\n");

  free_2d_array(tracer_result);
  brk5_finalize(&brk45_data);
  free(prev_point);
  free(point);
  free(xp);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);  
}


double estimate_R_from_z_psi(double psi_value, double r0, double z0,  Equilibrium *equ)
{
  double delta_r = 0.2;
  double left_psi, right_psi;
  double tol = 1.0E-12;
  int max_iter = 500;
  double x0 = r0 - delta_r;
  double x1 = r0 + delta_r;

  cubicherm2d1f(x0, z0, equ->nw, equ->r, equ->nh, equ->z,
                equ->psi, &(left_psi), NULL, NULL, NULL);
  cubicherm2d1f(x1, z0, equ->nw, equ->r, equ->nh, equ->z,
                equ->psi, &(right_psi), NULL, NULL, NULL);

  if (psi_value >= max(left_psi, right_psi) || psi_value <= min(left_psi, right_psi) )
  {
    printf("the psi value is out of range!\n");
    exit(1);
  }
// TEST: use secant method to fine the point (r,z0) that the psi is psi_value
  double f0 = left_psi - psi_value;
  double f1 = right_psi - psi_value;

  if (fabs(f0)  < tol)
  {
    return x0;
  }
  if (fabs(f1) < tol)
  {
    return x1;
  }

  for (int iter = 0; iter < max_iter; iter++) 
  {
    if (fabs(f1-f0) < tol) 
    {
      printf("Error: Division by zero in iteration %d\n", iter);
      return NAN;
    }
    double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x0 = x1;
    f0 = f1;
    x1 = x2;
    cubicherm2d1f(x1, z0, equ->nw, equ->r, equ->nh, equ->z,
                equ->psi, &(f1), NULL, NULL, NULL);
    f1 = f1-psi_value;
    
    if (fabs(f1) < tol)
    {
      printf("Secant iteration: %d\n", iter);
      return x1;
    }
  }
  printf("Warning: Secant method did not converge after %d iterations.\n", max_iter);
  return NAN;
}
