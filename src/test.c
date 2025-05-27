#include "test.h"
#include "w3lib.h"

// Test module
void separatrix_test()
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

  free(xp);
  free_equilibrium(&dtt_example);  
}

static void f(double* x, double* fx)
{
  double tmp = *x;
  *fx=1.0*tmp*tmp*tmp + 2.5*tmp*tmp + 5.0*tmp +6.0;
}
static void f_dfdx(double* x, double* dfdx)
{
  double tmp = *x;
  *dfdx=3.0*tmp*tmp + 5.0*tmp + 5.0;
}

void interp1d_test()
{
  int n = 20;
  double x[n];
  double fx[n];
  double dfdx[n];
  for(int i=0; i < n; i++)
  {
    x[i] = i;
    f(&x[i], &fx[i]);
    f_dfdx(&x[i], &dfdx[i]);
  }
  Interp1DFunction* interp;
  interp = create_cubicherm1D_interp(x,fx,dfdx,n);
  double dx=0.2;
  int m = (int)(1.0 / dx) * (n-1);
  printf("m: %d\n",m);
  double x1[m];
  double y1[m];
  double y2[m];

  for(int i=0; i < m; i++)
  {
    x1[i]=dx*(double)(i);
    interp->eval(interp->data,x1[i],&y1[i]);
    printf("%.15f %.15f \n", x1[i], y1[i]);
    f(&x1[i], &y2[i]);

  }
  const char *filename= "interp1d_test";
  FILE* file = fopen(filename, "w");
  for (int i = 0; i<m; i++)
  {
    // printf("%.15f %.15f %.15f\n", x1[i], y2[i], y1[i]);
    fprintf(file, "%.15f %.15f %.15f\n", x1[i], y2[i], y1[i]);
  }
  fclose(file);
  printf("Finish interp_test.\n");
  free_interp1d_function(interp);
}

void test_find_xpoint()
{
  printf("*******************************************\n");
  printf("Begint to test find_xpoint\n"); 

  InputPara w3_input;
  init_inputpara(&w3_input);
  print_inputpara(&w3_input);
  Equilibrium dtt_example;
  
  init_equilibrium(&dtt_example);

  read_equilib_geqdsk(&dtt_example,w3_input.equilibrium_file);
  print_equilibrium(&dtt_example);
  
  int xpt_n = 2;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = 1.85;
  est_xpt[0][1] = -1.16;

  est_xpt[1][0] = 1.58;
  est_xpt[1][1] = 1.61;

  interpl_1D_f interpl_1D_f = cubicherm2d1f;
  interpl_2D_f interpl_2D_f = cubicherm2d2f;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_1D_f, interpl_2D_f, xpt_array);

  free_2d_array(est_xpt);
  free_equilibrium(&dtt_example);  

}

void new_separatrix_test(){
  InputPara w3_input;
  init_inputpara(&w3_input);
  print_inputpara(&w3_input);
  Equilibrium dtt_example;
  
  init_equilibrium(&dtt_example);

  read_equilib_geqdsk(&dtt_example,w3_input.equilibrium_file);
  print_equilibrium(&dtt_example);
  
  int xpt_n = 2;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = 1.85;
  est_xpt[0][1] = -1.16;

  est_xpt[1][0] = 1.58;
  est_xpt[1][1] = 1.61;

  interpl_1D_f interpl_1D_f = cubicherm2d1f;
  interpl_2D_f interpl_2D_f = cubicherm2d2f;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_1D_f, interpl_2D_f, xpt_array);

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);

  // write_mag_field_torsys(&test_magfield);


  SeparatrixStr* sep=init_separatrix_default();
  void *p1;
  void *p2;
  p1=NULL;
  p2=NULL;

//build the interpolator; x_tmp,fx_tmp, dfdx_tmp are nothing realted to x or y. 

  Interp1DFunction* interp=create_cubicherm1D_interp(NULL, NULL, NULL, 2);

/************************************************
*  Build the tracer for generation separatrix   *
************************************************/ 
  double direction[3]={1.0,1.0,1.0};
  RKSolverData brk45_data;

  double stepsize = 0.1;

  ode_function ode_func = {
    .ndim = 2,
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
  brk45_solver.initialize(&brk45_data);

  generate_separatrix_bytracing(sep, &xpt_array[1], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);

  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_separatrix_default(sep);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);

};