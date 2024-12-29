/* This code is just test the magnetic field line tracer.
 * Should not be included in source code. 
*/
#include "tracertest.h"

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
  
  const char *filename4="tracer_rror";
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
  
  
  
  free(x1);
  free(y_euler);
  free(y_brk5);
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

  int nx = 101;
  int ny = 201;
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
  int nx1 = 301;
  int ny1 = 601;
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
      cubicherm_2d(x1[i],y1[j],nx, x, ny, y, f, f0, f1, NULL, NULL, NULL);
      // printf("finish cubicherm_2d interpolator\n");
      f0 = &f_bilinear[i][j][0];
      f1 = &f_bilinear[i][j][1];
      bilenar_2d(x1[i],y1[j],nx, x, ny, y, f, f0, f1, NULL, NULL, NULL);
      // printf("finish bilenar_2d interpolator\n");
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
      fprintf(file1, "%.15f %.15f %.15f\n", x[i], y[j], f[i][j][1]);
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
  double x = 1.534;
  double y = -1.125;
  value = get_psi_from_rz(&dtt_example,x, y);
  printf("psi at %lf %lf is %lf\n",x, y, value);


//**********Test four lines of the separatrix tracing*******************/

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
//**********************************************************************/


/*test magnetic field line calculation*/

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_2nd";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);
  write_mag_field_torsys(&test_magfield);
  free(xp);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);  
}
