#include "test.h"
#include "w3lib.h"
#include <time.h>

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
  xp = find_Xpoint_from_DG(&dtt_example, w3_input.xpt_estimation);

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
  write_DLList(line_list[0],"sep0");
  write_DLList(line_list[1],"sep1");
  write_DLList(line_list[2],"sep2");
  write_DLList(line_list[3],"sep3");

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

  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

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

  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);

  // write_mag_field_torsys(&test_magfield);




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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[1], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);
  
  free_separatrix_default(sep);


  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);

};

void read_trg_test()
{
  char* name="example.trg";
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, name);
  free_dgtrg(trg);
}

void get_target_curve_from_trg_test()
{
  char* name="example.trg";
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, name);
  int n=0;
  TargetDLListCurve* tgt_cur1=create_target_curve_from_dgtrg(trg,n);
  printf("%d\n" ,trg->n_target_curve[0]);
  // printf("%d\n" ,trg->n_target_curve[1]);
  // // TargetDLListCurve* tgt_cur2=create_target_curve_from_dgtrg(trg,2);
  printf_target_curve(tgt_cur1);
  // // printf_target_curve(tgt_cur2);
  free_target_curve(tgt_cur1);
  // free_target_curve(tgt_cur2);
  free_dgtrg(trg);

}
void target_curve_test()
{
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

  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);

  // write_mag_field_torsys(&test_magfield);




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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[1], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);



  char* name="example.trg";
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, name);
  int n=0;
  TargetDLListCurve* tgt_cur=create_target_curve_from_dgtrg(trg,n);
  // printf("%d\n" ,trg->n_target_curve[0]);
  // printf("%d\n" ,trg->n_target_curve[1]);
  // // TargetDLListCurve* tgt_cur2=create_target_curve_from_dgtrg(trg,2);
  // printf_target_curve(tgt_cur);
  // // printf_target_curve(tgt_cur2);
  // free_target_curve(tgt_cur1);
  // free_target_curve(tgt_cur2);
  
  DLListNode* tgt_cur_head=tgt_cur->head;
  double itrsct_r, itrsct_z;
  int num_itrsct;
  for(int i=0; i<4; i++)
  {
    int value;
    value = has_intersection_DLList(tgt_cur_head, sep->line_list[i]);
    if(value==0)
    {
      printf("The %d sep line intersects with target cureve.\n",i);
      insert_intersections_DLList(tgt_cur_head, sep->line_list[i], &itrsct_r, &itrsct_z);
      printf("intersection:\n");
      printf("%lf %lf\n", itrsct_r,itrsct_z);
      num_itrsct=i;
    }
    else
    {
      //printf("The %d sep line do NOT intersect with target cureve.\n", i);
    }
  }
  int n_cut;
  n_cut = cut_DLList_after_point(sep->line_list[num_itrsct], itrsct_r,itrsct_z);
  printf("Deleted points: %d\n", n_cut);
  char* cut_sep="sep_cut_baseline";
  write_DLList(sep->line_list[num_itrsct], cut_sep);
  char* name_tgt_cur="target_curve";
  write_DLList(tgt_cur_head, name_tgt_cur);




  
  free_target_curve(tgt_cur);
  free_dgtrg(trg);
  free_separatrix_default(sep);
  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);
}


void grad_psi_test()
{

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

  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);


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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[1], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);
  


  GradPsiStr *gradpsi=init_grad_psi();
  calc_grad_psi(&dtt_example, gradpsi, central_diff_2nd_2d);
  char name[32]="gradpsi";
  write_grad_psi(gradpsi, name);


// change the odf function for gradpsi line tracing
  ode_func.compute_f=ode_f_gradpsi_cubicherm;
  ode_func.data=gradpsi;

  GradPsiLineStr* gradpsilines=init_gradpsiline_default();
  
  //Create a opoint structure
  OPointStr* opoint = create_opoint();
  opoint->centerX=2.27;
  opoint->centerY=0.187;
  
  generate_gradpsiline_bytracing(gradpsilines, gradpsi, opoint, sep, NULL, &ode_func, &brk45_solver);


  free_gradpsiline_default(gradpsilines);
  free_grad_psi(gradpsi);
  free_separatrix_default(sep);

  free_opoint(opoint);
  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);
}


void divgeo_test()
{
  W3Config w3config;
  load_w3_config("w3.ini",&w3config);
  print_w3_config(&w3config);
  
  Equilibrium dtt_example;
  init_equilibrium(&dtt_example);
  read_equilib_geqdsk(&dtt_example,w3config.file_config.geqdsk_file);
  print_equilibrium(&dtt_example);


  int xpt_n = w3config.grid2d_config.xpoint_number;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = w3config.grid2d_config.xpoint_r_est[xpt_n-1];
  est_xpt[0][1] = w3config.grid2d_config.xpoint_z_est[xpt_n-1];


  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[1];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);


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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[0], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);
  
  GradPsiStr *gradpsi=init_grad_psi();
  calc_grad_psi(&dtt_example, gradpsi, central_diff_2nd_2d);
  char name[32]="gradpsi";
  write_grad_psi(gradpsi, name);

// change the odf function for gradpsi line tracing
  ode_func.compute_f=ode_f_gradpsi_cubicherm;
  ode_func.data=gradpsi;

  GradPsiLineStr* gradpsilines=init_gradpsiline_default();
  
  //Create a opoint structure
  OPointStr* opoint = create_opoint();
  opoint->centerX=2.27;
  opoint->centerY=0.18;
  
  generate_gradpsiline_bytracing(gradpsilines, gradpsi, opoint, sep, NULL, &ode_func, &brk45_solver);


  char* trgname=w3config.file_config.divgeo_trg_file;
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, trgname);

/***********************************************
*   Update trg regions(psi valuse)
***********************************************/
  for(int i=0; i<3; i++)
  {
    trg->regions[i]->level[0]=xpt_array[1].level;
  }

  write_sn_gridzoneinfo_from_dgtrg(trg, &dtt_example, sep, gradpsilines, &w3config.grid2d_config);



/***********************************************
*   Read input of GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_SOL");
  GridZoneInfo* pfrgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_PFR");
  GridZoneInfo* coregzinfo=load_GridZoneInfo_from_input("gridzoneinfo_CORE");

  print_GridZoneInfo(solgzinfo);
  print_GridZoneInfo(pfrgzinfo);
  print_GridZoneInfo(coregzinfo);


  write_polsegms_from_dgtrg(trg,"polseginfo");
  PolSegmsInfo* polseginfo=read_PolSegmsInfo_from_file("polseginfo");
  print_PolSegmsInfo(polseginfo);
  free_PolSegmsInfo(polseginfo);
  
  free_GridZoneInfo(&solgzinfo);
  free_GridZoneInfo(&pfrgzinfo);
  free_GridZoneInfo(&coregzinfo);

  free_dgtrg(trg);
  free_gradpsiline_default(gradpsilines);
  free_grad_psi(gradpsi);
  free_separatrix_default(sep);

  free_opoint(opoint);
  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);
}


void meshgeneration_test()
{
  InputPara w3_input;
  init_inputpara(&w3_input);
  print_inputpara(&w3_input);
  Equilibrium dtt_example;
  init_equilibrium(&dtt_example);
  read_equilib_geqdsk(&dtt_example,w3_input.equilibrium_file);
  correct_direction_lower_divertor(&dtt_example);
  print_equilibrium(&dtt_example);


  int xpt_n = 2;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = 1.85;
  est_xpt[0][1] = -1.16;

  est_xpt[1][0] = 1.58;
  est_xpt[1][1] = 1.61;

  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);


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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[1], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);
  
  GradPsiStr *gradpsi=init_grad_psi();
  calc_grad_psi(&dtt_example, gradpsi, central_diff_2nd_2d);
  char name[32]="gradpsi";
  write_grad_psi(gradpsi, name);

// change the odf function for gradpsi line tracing
  ode_func.compute_f=ode_f_gradpsi_cubicherm;
  ode_func.data=gradpsi;

  GradPsiLineStr* gradpsilines=init_gradpsiline_default();
  
  //Create a opoint structure
  OPointStr* opoint = create_opoint();
  opoint->centerX=2.27;
  opoint->centerY=0.187;
  
  generate_gradpsiline_bytracing(gradpsilines, gradpsi, opoint, sep, NULL, &ode_func, &brk45_solver);

/***********************************************
*   1. Read input of GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_SOL");
  GridZoneInfo* pfrgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_PFR");
  GridZoneInfo* coregzinfo=load_GridZoneInfo_from_input("gridzoneinfo_CORE");

/***********************************************
*   2. Read input of polsegminfo
***********************************************/
  PolSegmsInfo* polseginfo=read_PolSegmsInfo_from_file("polseginfo");

  
/***********************************************
*   3. Create separatrix distribution
***********************************************/
  SepDistStr* sepdist=create_SepDistStr_from_sep(sep);
  update_sn_SepDistStr_from_GridZoneInfo(sepdist,solgzinfo);
  update_sn_SepDistStr_from_PolSegmsInfo(sepdist,polseginfo);
  update_SepDistStr_gridpoint_curve(sepdist);
  int idx=sepdist->index[0];
  write_curve("gridpoint_curve0", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[1];
  write_curve("gridpoint_curve1", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[2];
  write_curve("gridpoint_curve2", sepdist->edges[idx]->gridpoint_curve);

  // write_DLList(sepdist->edges[sepdist->index[0]]->head,"sepdist_line0");
  // write_DLList(sepdist->edges[sepdist->index[1]]->head,"sepdist_line1");
  // write_DLList(sepdist->edges[sepdist->index[2]]->head,"sepdist_line2");
  // write_DLList(sepdist->edges[sepdist->index[3]]->head,"sepdist_line3");

/***********************************************
*   4. Create gridzone
***********************************************/
  GridZone* solgz=create_sn_CARRE2D_GridZone(solgzinfo, sepdist);
  GridZone* pfrgz=create_sn_CARRE2D_GridZone(pfrgzinfo, sepdist);
  GridZone* coregz=create_sn_CARRE2D_GridZone(coregzinfo, sepdist);

  write_curve("sol_gz_c",solgz->first_bnd_curve);
  write_curve("sol_gz_gpc",solgz->first_gridpoint_curve);
  write_curve("pfr_gz_c",pfrgz->first_bnd_curve);
  write_curve("pfr_gz_gpc",pfrgz->first_gridpoint_curve);
  write_curve("core_gz_c",coregz->first_bnd_curve);
  write_curve("core_gz_gpc",coregz->first_gridpoint_curve);

  
  
/***********************************************
*   5. Generater grid for each gridzone
***********************************************/
  
  // change the odf function for gradpsi line tracing
  ode_func.compute_f=ode_f_brz_torsys_cubicherm;
  ode_func.data=&test_magfield;

  
  TwoDimGrid* sol2dgrid=create_2Dgrid_poloidal_major(solgz->first_gridpoint_curve->n_point, solgz->nr);
  generate_CARRE_2Dgrid_default(sol2dgrid, solgz, &ode_func, &brk45_solver);
  generate_CARRE_2Dgrid_default(sol2dgrid, pfrgz, &ode_func, &brk45_solver);
  generate_CARRE_2Dgrid_default(sol2dgrid, coregz, &ode_func, &brk45_solver);

  


/***********************************************
*   6. Free space
***********************************************/

  free_2Dgrid(sol2dgrid);
  
  free_GridZone(solgz);
  free_GridZone(pfrgz);
  free_GridZone(coregz);

  free_GridZoneInfo(&solgzinfo);
  free_GridZoneInfo(&pfrgzinfo);
  free_GridZoneInfo(&coregzinfo);

  free_PolSegmsInfo(polseginfo);

  free_SepDistStr(sepdist);
  free_gradpsiline_default(gradpsilines);
  free_grad_psi(gradpsi);
  free_separatrix_default(sep);

  free_opoint(opoint);
  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);
}


void ThreeDimMeshGeneration_test()
{

  InputPara w3_input;
  init_inputpara(&w3_input);
  print_inputpara(&w3_input);
  Equilibrium dtt_example;
  init_equilibrium(&dtt_example);
  read_equilib_geqdsk(&dtt_example,w3_input.equilibrium_file);
  correct_direction_lower_divertor(&dtt_example);
  print_equilibrium(&dtt_example);


  int xpt_n = 1;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = 1.86;
  est_xpt[0][1] = -1.16;


  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[1];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);


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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[0], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);
  
  GradPsiStr *gradpsi=init_grad_psi();
  calc_grad_psi(&dtt_example, gradpsi, central_diff_2nd_2d);
  char name[32]="gradpsi";
  write_grad_psi(gradpsi, name);

// change the odf function for gradpsi line tracing
  ode_func.compute_f=ode_f_gradpsi_cubicherm;
  ode_func.data=gradpsi;

  GradPsiLineStr* gradpsilines=init_gradpsiline_default();
  
  //Create a opoint structure
  OPointStr* opoint = create_opoint();
  opoint->centerX=2.27;
  opoint->centerY=0.18;
  
  generate_gradpsiline_bytracing(gradpsilines, gradpsi, opoint, sep, NULL, &ode_func, &brk45_solver);

/***********************************************
*   0. Create GridZoneInfo and polseginfo
***********************************************/

  // char* trgname="example.trg";
  // DivGeoTrg* trg=create_dgtrg();
  // int status=load_dgtrg_from_file(trg, trgname);

  // //because the is numeric error of psi in X-point.
  // //The psi of strike point may not identical with X-point.
  // //Make them the same
  // for(int i=0; i<3; i++)
  // {
  //   trg->regions[i]->level[0]=xpt_array[1].level;
  // }

  // write_sn_gridzoneinfo_from_dgtrg(trg, &dtt_example, sep, gradpsilines);

/***********************************************
*   1. Read input of GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_SOL");
  GridZoneInfo* pfrgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_PFR");
  GridZoneInfo* coregzinfo=load_GridZoneInfo_from_input("gridzoneinfo_CORE");

  print_GridZoneInfo(solgzinfo);
  print_GridZoneInfo(pfrgzinfo);
  print_GridZoneInfo(coregzinfo);

/***********************************************
*   2. Read input of polsegminfo
***********************************************/
  PolSegmsInfo* polseginfo=read_PolSegmsInfo_from_file("polseginfo");

  
/***********************************************
*   3. Create separatrix distribution
***********************************************/
  SepDistStr* sepdist=create_SepDistStr_from_sep(sep);
  update_sn_SepDistStr_from_GridZoneInfo(sepdist,solgzinfo);
  update_sn_SepDistStr_from_PolSegmsInfo(sepdist,polseginfo);
  update_SepDistStr_gridpoint_curve(sepdist);
  int idx=sepdist->index[0];
  write_curve("gridpoint_curve0", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[1];
  write_curve("gridpoint_curve1", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[2];
  write_curve("gridpoint_curve2", sepdist->edges[idx]->gridpoint_curve);

  // write_DLList(sepdist->edges[sepdist->index[0]]->head,"sepdist_line0");
  // write_DLList(sepdist->edges[sepdist->index[1]]->head,"sepdist_line1");
  // write_DLList(sepdist->edges[sepdist->index[2]]->head,"sepdist_line2");
  // write_DLList(sepdist->edges[sepdist->index[3]]->head,"sepdist_line3");


/***********************************************
*   4. Update polseparatrix distribution 
***********************************************/
  ode_func.compute_f=ode_f_brz_torsys_cubicherm;
  ode_func.data=&test_magfield;
  ode_func.ndim=3;

  /*=================================
   Toroidal phi defintion
  ==================================*/
  int nphi=11;
  double delta=2;
  const int idx_mid=nphi/2;


  double *phi=malloc(nphi*sizeof(double));
  phi[0]=10.0;
  for(int i=1;i<nphi;i++)
  {
    phi[i]=phi[i-1]+delta;
  }
  // phi[4]=30;
  // phi[5]=32;
  // phi[6]=34;

  update_sn_SepDistStr_PolSegmsInfo_EMC3_2Dgrid(polseginfo, sepdist, &ode_func, &brk45_solver,
                                                phi[idx_mid], nphi, phi);
  write_PolSegmsInfo(polseginfo, "polseginfo_3DGRID");

/***********************************************
*   5. Create gridzone and 2D basegrid
***********************************************/
  GridZone* solgz=create_sn_CARRE2D_GridZone(solgzinfo, sepdist);
  GridZone* pfrgz=create_sn_CARRE2D_GridZone(pfrgzinfo, sepdist);
  GridZone* coregz=create_sn_CARRE2D_GridZone(coregzinfo, sepdist);

  TwoDimGrid* sol2dgrid=create_2Dgrid_poloidal_major(solgz->first_gridpoint_curve->n_point, solgz->nr);
  generate_EMC3_2Dgrid_default(sol2dgrid, solgz, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);

  TwoDimGrid* pfr2dgrid=create_2Dgrid_poloidal_major(pfrgz->first_gridpoint_curve->n_point, pfrgz->nr);
  generate_EMC3_2Dgrid_default(pfr2dgrid, pfrgz, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);

  //Becaure full about CORE GRIDZONE
  TwoDimGrid* core2dgrid=create_2Dgrid_poloidal_major(coregz->first_gridpoint_curve->n_point, coregz->nr);
  ode_func.ndim=2;
  generate_CARRE_2Dgrid_default(core2dgrid, coregz, &ode_func, &brk45_solver);
  ode_func.ndim=3;

//   write_curve("sol_gz_c",solgz->first_bnd_curve);
//   write_curve("sol_gz_gpc",solgz->first_gridpoint_curve);
//   write_curve("pfr_gz_c",pfrgz->first_bnd_curve);
//   write_curve("pfr_gz_gpc",pfrgz->first_gridpoint_curve);
//   write_curve("core_gz_c",coregz->first_bnd_curve);
//   write_curve("core_gz_gpc",coregz->first_gridpoint_curve);

  
  
/**********************************************
*   6. Expand the 2D basegrid                 *
***********************************************/
  TwoDimGrid* sol2dgrid_exptgt=load_2Dgrid_from_file("SOLGRIDZONEINFO_3D_2DBASEGRID");
  TwoDimGrid* pfr2dgrid_exptgt=load_2Dgrid_from_file("PFRGRIDZONEINFO_3D_2DBASEGRID");

  expand_target_EMC3_2Dgrid_default(sol2dgrid_exptgt, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);
  write_2Dgrid(sol2dgrid_exptgt,"SOLGRIDZONEINFO_3D_2DEXPTGT");
 
  expand_target_EMC3_2Dgrid_default(pfr2dgrid_exptgt, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);
  write_2Dgrid(pfr2dgrid_exptgt,"PFRGRIDZONEINFO_3D_2DEXPTGT");

/**********************************************
*   7. generate EMC3 3D GRID                  *
***********************************************/
  TwoDimGrid* grid_tmp1=create_2Dgrid_poloidal_major(sol2dgrid_exptgt->npol, sol2dgrid_exptgt->nrad);
  TwoDimGrid* grid_tmp2=create_2Dgrid_poloidal_major(pfr2dgrid_exptgt->npol, pfr2dgrid_exptgt->nrad);
  TwoDimGrid* grid_tmp3=create_2Dgrid_poloidal_major(core2dgrid->npol, core2dgrid->nrad);


  // restore_3D_mag_direction(&ode_func);
  // for(int i=0;i<1;i++)
  // {
  //   char name_tmp[32];
  //   int m = 1+(i+1);
  //   sprintf(name_tmp,"SOL_3D_%i",m);
  //   generate_2Dgrid_tracing(sol2dgrid_exptgt, phi[idx_mid], grid_tmp1, phi[m], &ode_func, &brk45_solver);
  //   write_2Dgrid_RZCSYS_to_XYZCSYS(grid_tmp1, phi[m], name_tmp);
  // }

  // char name_tmp[32];
  // sprintf(name_tmp,"SOL_3D_1"); 
  // write_2Dgrid_RZCSYS_to_XYZCSYS(sol2dgrid_exptgt, phi[idx_mid], name_tmp);

  // restore_3D_mag_direction(&ode_func);
  // for(int i=0;i<1;i++)
  // {
  //   char name_tmp[32];
  //   int m = 1-(i+1);
  //   sprintf(name_tmp,"SOL_3D_%i",m);
  //   generate_2Dgrid_tracing(sol2dgrid_exptgt, phi[idx_mid], grid_tmp1, phi[m], &ode_func, &brk45_solver);
  //   write_2Dgrid_RZCSYS_to_XYZCSYS(grid_tmp1, phi[m], name_tmp);
  // }


  restore_3D_mag_direction(&ode_func);
  ThreeDimGrid* sol3dgrid=create_3Dgrid_radial_major(sol2dgrid_exptgt->npol,sol2dgrid_exptgt->nrad,nphi);
  ThreeDimGrid* pfr3dgrid=create_3Dgrid_radial_major(pfr2dgrid_exptgt->npol,pfr2dgrid_exptgt->nrad,nphi);
  ThreeDimGrid* core3dgrid=create_3Dgrid_radial_major(core2dgrid->npol,core2dgrid->nrad,nphi);

  generate_EMC3_3Dgrid_from_2Dgrid_tracing(sol2dgrid_exptgt, sol3dgrid, phi[idx_mid], nphi, phi, &ode_func, &brk45_solver);
  generate_EMC3_3Dgrid_from_2Dgrid_tracing(pfr2dgrid_exptgt, pfr3dgrid, phi[idx_mid], nphi, phi, &ode_func, &brk45_solver);
  generate_EMC3_3Dgrid_from_2Dgrid_tracing(core2dgrid, core3dgrid, phi[idx_mid], nphi, phi, &ode_func, &brk45_solver);


  write_EMC3_3Dgrid_to_XYZ_CSYS(sol3dgrid,"EMC3_3DGRID_SOL");
  write_EMC3_3Dgrid_to_EMC3_format(sol3dgrid,"grid3D_SOL.dat",false,false);

  write_EMC3_3Dgrid_to_XYZ_CSYS(core3dgrid,"EMC3_3DGRID_CORE");
  write_EMC3_3Dgrid_to_EMC3_format(core3dgrid,"grid3D_CORE.dat",false,false);

  write_EMC3_3Dgrid_to_XYZ_CSYS(pfr3dgrid,"EMC3_3DGRID_PFR");
  write_EMC3_3Dgrid_to_EMC3_format(pfr3dgrid,"grid3D_PFR.dat",false,false);
  
/**********************************************
*   8. Free space                             *
***********************************************/
  free(phi);
  free_3Dgrid(sol3dgrid);
  free_2Dgrid(grid_tmp1);
  free_2Dgrid(grid_tmp2);
  free_2Dgrid(grid_tmp3);

  free_2Dgrid(sol2dgrid_exptgt);
  free_2Dgrid(pfr2dgrid_exptgt);


  free_2Dgrid(sol2dgrid);
  free_2Dgrid(core2dgrid);

  
  free_GridZone(solgz);
  free_GridZone(pfrgz);
  free_GridZone(coregz);

  free_GridZoneInfo(&solgzinfo);
  free_GridZoneInfo(&pfrgzinfo);
  free_GridZoneInfo(&coregzinfo);

  free_PolSegmsInfo(polseginfo);


  free_SepDistStr(sepdist);
  free_gradpsiline_default(gradpsilines);
  free_grad_psi(gradpsi);
  free_separatrix_default(sep);

  free_opoint(opoint);
  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);
}



void EMC3_3D_grid_generation_test()
{
/**************************
*  Read W3 Configuration  *
**************************/ 
  W3Config w3config;
  load_w3_config("w3.ini",&w3config);
  print_w3_config(&w3config);
  // ConfigData raw_config = {0};
  // parse_config_file("w3.ini",&raw_config);
  // print_config_data(&raw_config);
  printf("Press Enter to continue...");
  getchar();
  printf("Program continues...\n");
  
  Equilibrium dtt_example;
  init_equilibrium(&dtt_example);
  read_equilib_geqdsk(&dtt_example,w3config.file_config.geqdsk_file);
  correct_direction_lower_divertor(&dtt_example);
  print_equilibrium(&dtt_example);
  int xpt_n = w3config.grid2d_config.xpoint_number;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = w3config.grid2d_config.xpoint_r_est[xpt_n-1];
  est_xpt[0][1] = w3config.grid2d_config.xpoint_z_est[xpt_n-1];

  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[1];
  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  
  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);
  write_mag_field_torsys(&test_magfield);

  printf("Press Enter to continue...");
  getchar();
  printf("Program continues...\n");

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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[0], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);
  
  GradPsiStr *gradpsi=init_grad_psi();
  calc_grad_psi(&dtt_example, gradpsi, central_diff_2nd_2d);
  char name[32]="gradpsi";
  write_grad_psi(gradpsi, name);

  printf("Press Enter to continue...");
  getchar();
  printf("Program continues...\n");

// change the odf function for gradpsi line tracing
  ode_func.compute_f=ode_f_gradpsi_cubicherm;
  ode_func.data=gradpsi;

  GradPsiLineStr* gradpsilines=init_gradpsiline_default();
  
  //Create a opoint structure
  OPointStr* opoint = create_opoint();
  opoint->centerX=2.27;
  opoint->centerY=0.18;
  
  generate_gradpsiline_bytracing(gradpsilines, gradpsi, opoint, sep, NULL, &ode_func, &brk45_solver);


  char* trgname=w3config.file_config.divgeo_trg_file;
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, trgname);

/***********************************************
*   Update trg regions(psi valuse)
***********************************************/
// ensure the first point is exact same with X-point.
  for(int i=0; i<3; i++)
  {
    trg->regions[i]->level[0]=xpt_array[1].level;
  }

// !! Correct psi values according to the sibry
  correct_psi_based_on_sibry(&dtt_example);

  write_sn_gridzoneinfo_from_dgtrg(trg, &dtt_example, sep, gradpsilines, &w3config.grid2d_config);

  write_polsegms_from_dgtrg(trg,"polseginfo");


/***********************************************
*   1. Read input of GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_SOL");
  GridZoneInfo* pfrgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_PFR");
  GridZoneInfo* coregzinfo=load_GridZoneInfo_from_input("gridzoneinfo_CORE");

  print_GridZoneInfo(solgzinfo);
  print_GridZoneInfo(pfrgzinfo);
  print_GridZoneInfo(coregzinfo);

/***********************************************
*   2. Read input of polsegminfo
***********************************************/
  PolSegmsInfo* polseginfo=read_PolSegmsInfo_from_file("polseginfo");




/***********************************************
*   3. Create separatrix distribution
***********************************************/
  SepDistStr* sepdist=create_SepDistStr_from_sep(sep);
  update_sn_SepDistStr_from_GridZoneInfo(sepdist,solgzinfo);
  update_sn_SepDistStr_from_PolSegmsInfo(sepdist,polseginfo);
  update_SepDistStr_gridpoint_curve(sepdist);
  int idx=sepdist->index[0];
  write_curve("gridpoint_curve0", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[1];
  write_curve("gridpoint_curve1", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[2];
  write_curve("gridpoint_curve2", sepdist->edges[idx]->gridpoint_curve);



/***********************************************
*   4. Update polseparatrix distribution 
***********************************************/
  ode_func.compute_f=ode_f_brz_torsys_cubicherm;
  ode_func.data=&test_magfield;
  ode_func.ndim=3;

  /*=================================
   Toroidal phi defintion
  ==================================*/
  //points number = cells number + 1
  int nphi=w3config.grid3d_config.toroidal_cell_number+1;
  double delta=w3config.grid3d_config.toroidal_delta;
  const int idx_mid=nphi/2;


  double *phi=malloc(nphi*sizeof(double));
  phi[0]=0.0;
  for(int i=1;i<nphi;i++)
  {
    phi[i]=phi[i-1]+delta;
  }
  // phi[4]=30;
  // phi[5]=32;
  // phi[6]=34;

  update_sn_SepDistStr_PolSegmsInfo_EMC3_2Dgrid(polseginfo, sepdist, &ode_func, &brk45_solver,
                                                phi[idx_mid], nphi, phi);
  write_PolSegmsInfo(polseginfo, "polseginfo_3DGRID");


/***********************************************
*   5. Create gridzone and 2D basegrid
***********************************************/
  GridZone* solgz=create_sn_CARRE2D_GridZone(solgzinfo, sepdist);
  GridZone* pfrgz=create_sn_CARRE2D_GridZone(pfrgzinfo, sepdist);
  GridZone* coregz=create_sn_CARRE2D_GridZone(coregzinfo, sepdist);

  TwoDimGrid* sol2dgrid=create_2Dgrid_poloidal_major(solgz->first_gridpoint_curve->n_point, solgz->nr);
  generate_EMC3_2Dgrid_default(sol2dgrid, solgz, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);

  TwoDimGrid* pfr2dgrid=create_2Dgrid_poloidal_major(pfrgz->first_gridpoint_curve->n_point, pfrgz->nr);
  generate_EMC3_2Dgrid_default(pfr2dgrid, pfrgz, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);

  //Becaure full about CORE GRIDZONE
  TwoDimGrid* core2dgrid=create_2Dgrid_poloidal_major(coregz->first_gridpoint_curve->n_point, coregz->nr);
  ode_func.ndim=2;
  generate_CARRE_2Dgrid_default(core2dgrid, coregz, &ode_func, &brk45_solver);
  close_pol_first_last_2Dgrid(core2dgrid);
  ode_func.ndim=3;

/**********************************************
*   6. Expand the 2D basegrid                 *
***********************************************/
  TwoDimGrid* sol2dgrid_exptgt=load_2Dgrid_from_file("SOLGRIDZONEINFO_3D_2DBASEGRID");
  TwoDimGrid* pfr2dgrid_exptgt=load_2Dgrid_from_file("PFRGRIDZONEINFO_3D_2DBASEGRID");

  expand_target_EMC3_2Dgrid_default(sol2dgrid_exptgt, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);
  write_2Dgrid(sol2dgrid_exptgt,"SOLGRIDZONEINFO_3D_2DEXPTGT");
 
  expand_target_EMC3_2Dgrid_default(pfr2dgrid_exptgt, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);
  write_2Dgrid(pfr2dgrid_exptgt,"PFRGRIDZONEINFO_3D_2DEXPTGT");

/**********************************************
*   7. generate EMC3 3D GRID                  *
***********************************************/
  TwoDimGrid* grid_tmp1=create_2Dgrid_poloidal_major(sol2dgrid_exptgt->npol, sol2dgrid_exptgt->nrad);
  TwoDimGrid* grid_tmp2=create_2Dgrid_poloidal_major(pfr2dgrid_exptgt->npol, pfr2dgrid_exptgt->nrad);
  TwoDimGrid* grid_tmp3=create_2Dgrid_poloidal_major(core2dgrid->npol, core2dgrid->nrad);

  restore_3D_mag_direction(&ode_func);
  ThreeDimGrid* sol3dgrid=create_3Dgrid_radial_major(sol2dgrid_exptgt->npol,sol2dgrid_exptgt->nrad,nphi);
  ThreeDimGrid* pfr3dgrid=create_3Dgrid_radial_major(pfr2dgrid_exptgt->npol,pfr2dgrid_exptgt->nrad,nphi);
  ThreeDimGrid* core3dgrid=create_3Dgrid_radial_major(core2dgrid->npol,core2dgrid->nrad,nphi);

  generate_EMC3_3Dgrid_from_2Dgrid_tracing(sol2dgrid_exptgt, sol3dgrid, phi[idx_mid], nphi, phi, &ode_func, &brk45_solver);
  generate_EMC3_3Dgrid_from_2Dgrid_tracing(pfr2dgrid_exptgt, pfr3dgrid, phi[idx_mid], nphi, phi, &ode_func, &brk45_solver);
  generate_EMC3_3Dgrid_from_2Dgrid_tracing(core2dgrid,      core3dgrid, phi[idx_mid], nphi, phi, &ode_func, &brk45_solver);


  write_EMC3_3Dgrid_to_XYZ_CSYS(sol3dgrid,"XYZ_3DGRID_SOL");
  write_EMC3_3Dgrid_to_EMC3_format(sol3dgrid,"grid3D_SOL.dat",false,false);

  write_EMC3_3Dgrid_to_XYZ_CSYS(core3dgrid,"XYZ_3DGRID_CORE");
  write_EMC3_3Dgrid_to_EMC3_format(core3dgrid,"grid3D_CORE.dat",false,false);

  write_EMC3_3Dgrid_to_XYZ_CSYS(pfr3dgrid,"XYZ_3DGRID_PFR");
  write_EMC3_3Dgrid_to_EMC3_format(pfr3dgrid,"grid3D_PFR.dat",false,false);


/**********************************************
*   8. Free space                             *
***********************************************/
  free(phi);
  free_3Dgrid(sol3dgrid);
  free_3Dgrid(pfr3dgrid);
  free_3Dgrid(core3dgrid);


  free_2Dgrid(grid_tmp1);
  free_2Dgrid(grid_tmp2);
  free_2Dgrid(grid_tmp3);

  free_2Dgrid(sol2dgrid_exptgt);
  free_2Dgrid(pfr2dgrid_exptgt);

  free_2Dgrid(sol2dgrid);
  free_2Dgrid(pfr2dgrid);
  free_2Dgrid(core2dgrid);

  free_GridZone(solgz);
  free_GridZone(pfrgz);
  free_GridZone(coregz);

  free_PolSegmsInfo(polseginfo);
  free_GridZoneInfo(&solgzinfo);
  free_GridZoneInfo(&pfrgzinfo);
  free_GridZoneInfo(&coregzinfo);

  free_dgtrg(trg);
  free_SepDistStr(sepdist);
  free_gradpsiline_default(gradpsilines);
  free_grad_psi(gradpsi);
  free_separatrix_default(sep);

  free_opoint(opoint);
  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);

}

void Expanded_2D_grid_generation_test()
{
/**************************
*  Read W3 Configuration  *
**************************/ 
  W3Config w3config;
  load_w3_config("w3.ini",&w3config);
  print_w3_config(&w3config);
  // ConfigData raw_config = {0};
  // parse_config_file("w3.ini",&raw_config);
  // print_config_data(&raw_config);
  printf("Press Enter to continue...");
  getchar();
  printf("Program continues...\n");
  
  Equilibrium dtt_example;
  init_equilibrium(&dtt_example);
  read_equilib_geqdsk(&dtt_example,w3config.file_config.geqdsk_file);
  correct_direction_lower_divertor(&dtt_example);
  print_equilibrium(&dtt_example);
  int xpt_n = w3config.grid2d_config.xpoint_number;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = w3config.grid2d_config.xpoint_r_est[xpt_n-1];
  est_xpt[0][1] = w3config.grid2d_config.xpoint_z_est[xpt_n-1];

  interpl_2D_1f interpl_2D_1f = cubicherm2d1f;
  interpl_2D_2f interpl_2D_2f = cubicherm2d2f;

  _XPointInfo xpt_array[1];
  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_2D_1f, interpl_2D_2f, xpt_array);

  printf("Press Enter to continue...");
  getchar();
  printf("Program continues...\n");

  
  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);


//build the interpolator; x_tmp,fx_tmp, dfdx_tmp are nothing realted to x or y. 

  Interp1DFunction* interp=create_cubicherm1D_interp(NULL, NULL, NULL, 2);

/************************************************
*  Build the tracer for generation separatrix   *
************************************************/ 
  double direction[3]={1.0,1.0,1.0};
  RKSolverData brk45_data;

  double stepsize = 0.5;

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

  SeparatrixStr* sep=init_separatrix_default();
  generate_separatrix_bytracing(sep, &xpt_array[0], &dtt_example,&test_magfield, interp,&ode_func, &brk45_solver);
  
  GradPsiStr *gradpsi=init_grad_psi();
  calc_grad_psi(&dtt_example, gradpsi, central_diff_2nd_2d);
  char name[32]="gradpsi";
  write_grad_psi(gradpsi, name);

// change the odf function for gradpsi line tracing
  ode_func.compute_f=ode_f_gradpsi_cubicherm;
  ode_func.data=gradpsi;

  GradPsiLineStr* gradpsilines=init_gradpsiline_default();
  
  //Create a opoint structure
  OPointStr* opoint = create_opoint();
  opoint->centerX=2.27;
  opoint->centerY=0.18;
  
  generate_gradpsiline_bytracing(gradpsilines, gradpsi, opoint, sep, NULL, &ode_func, &brk45_solver);


  char* trgname=w3config.file_config.divgeo_trg_file;
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, trgname);

/***********************************************
*   Update trg regions(psi valuse)
***********************************************/
  for(int i=0; i<3; i++)
  {
    trg->regions[i]->level[0]=xpt_array[1].level;
  }

  write_sn_gridzoneinfo_from_dgtrg(trg, &dtt_example, sep, gradpsilines, &w3config.grid2d_config);

  write_polsegms_from_dgtrg(trg,"polseginfo");


/***********************************************
*   1. Read input of GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_SOL");
  GridZoneInfo* pfrgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_PFR");
  GridZoneInfo* coregzinfo=load_GridZoneInfo_from_input("gridzoneinfo_CORE");

  print_GridZoneInfo(solgzinfo);
  print_GridZoneInfo(pfrgzinfo);
  print_GridZoneInfo(coregzinfo);

/***********************************************
*   2. Read input of polsegminfo
***********************************************/
  PolSegmsInfo* polseginfo=read_PolSegmsInfo_from_file("polseginfo");

/***********************************************
*   3. Create separatrix distribution
***********************************************/
  SepDistStr* sepdist=create_SepDistStr_from_sep(sep);
  update_sn_SepDistStr_from_GridZoneInfo(sepdist,solgzinfo);
  update_sn_SepDistStr_from_PolSegmsInfo(sepdist,polseginfo);
  update_SepDistStr_gridpoint_curve(sepdist);
  int idx=sepdist->index[0];
  write_curve("gridpoint_curve0", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[1];
  write_curve("gridpoint_curve1", sepdist->edges[idx]->gridpoint_curve);
  idx=sepdist->index[2];
  write_curve("gridpoint_curve2", sepdist->edges[idx]->gridpoint_curve);



/***********************************************
*   4. Update polseparatrix distribution 
***********************************************/
  ode_func.compute_f=ode_f_brz_torsys_cubicherm;
  ode_func.data=&test_magfield;
  ode_func.ndim=3;

  /*=================================
   Toroidal phi defintion
  ==================================*/
  int nphi=11;
  double delta=2;
  const int idx_mid=nphi/2;


  double *phi=malloc(nphi*sizeof(double));
  phi[0]=10.0;
  for(int i=1;i<nphi;i++)
  {
    phi[i]=phi[i-1]+delta;
  }
  // phi[4]=30;
  // phi[5]=32;
  // phi[6]=34;

  update_sn_SepDistStr_PolSegmsInfo_EMC3_2Dgrid(polseginfo, sepdist, &ode_func, &brk45_solver,
                                                phi[idx_mid], nphi, phi);
  write_PolSegmsInfo(polseginfo, "polseginfo_3DGRID");


/***********************************************
*   5. Create gridzone and 2D basegrid
***********************************************/
  GridZone* solgz=create_sn_CARRE2D_GridZone(solgzinfo, sepdist);
  GridZone* pfrgz=create_sn_CARRE2D_GridZone(pfrgzinfo, sepdist);

  TwoDimGrid* sol2dgrid=create_2Dgrid_poloidal_major(solgz->first_gridpoint_curve->n_point, solgz->nr);
  generate_EMC3_2Dgrid_default(sol2dgrid, solgz, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);

  TwoDimGrid* pfr2dgrid=create_2Dgrid_poloidal_major(pfrgz->first_gridpoint_curve->n_point, pfrgz->nr);
  generate_EMC3_2Dgrid_default(pfr2dgrid, pfrgz, &ode_func, &brk45_solver, phi[idx_mid], nphi, phi);

/***********************************************
*   6. Expand the SOL neutral grid
***********************************************/

  //6.1 build the four boundaries
  DLListNode* sol_neu_top_ddl = load_DLList_from_file("SOL_neu_bnd");
  
  double r_tmp = get_x_2Dgrid(sol2dgrid, 0, sol2dgrid->nrad-1);
  double z_tmp = get_y_2Dgrid(sol2dgrid, 0, sol2dgrid->nrad-1);
  DLListNode* sol_neu_bottom_ddl = create_DLListNode(r_tmp, z_tmp);
  DLListNode* tail_tmp = sol_neu_bottom_ddl;

  for(int i = 1; i<sol2dgrid->npol;i++)
  {
    double r_tmp = get_x_2Dgrid(sol2dgrid, i, sol2dgrid->nrad-1);
    double z_tmp = get_y_2Dgrid(sol2dgrid, i, sol2dgrid->nrad-1);
    add_DLListnode_at_tail(&tail_tmp, r_tmp, z_tmp);
  }

  DLListNode* sol_neu_left_ddl = load_DLList_from_file("inner_targetcurve");
  DLListNode* sol_neu_right_ddl = load_DLList_from_file("outer_targetcurve");

  r_tmp = get_x_2Dgrid(sol2dgrid,  sol2dgrid->npol-1,sol2dgrid->nrad-1);
  z_tmp = get_y_2Dgrid(sol2dgrid,  sol2dgrid->npol-1,sol2dgrid->nrad-1);

  if(insert_point_for_DLList(sol_neu_right_ddl, r_tmp, z_tmp))
  {
    fprintf(stderr, "Unexpected error: the point is not in the outer target.\n");
    exit(EXIT_FAILURE);
  }

  if(cut_DLList_before_point(&sol_neu_right_ddl, r_tmp, z_tmp)==0)
  {
    fprintf(stderr, "Unexpected error: cannot cut the point in the outer target.\n");
    exit(EXIT_FAILURE);
  }

  r_tmp = get_x_2Dgrid(sol2dgrid, 0, sol2dgrid->nrad-1);
  z_tmp = get_y_2Dgrid(sol2dgrid, 0, sol2dgrid->nrad-1);

  if(insert_point_for_DLList(sol_neu_left_ddl, r_tmp, z_tmp))
  {
    fprintf(stderr, "Unexpected error: the point is not in the inner target.\n");
    exit(EXIT_FAILURE);
  }

  if(cut_DLList_before_point(&sol_neu_left_ddl, r_tmp, z_tmp)==0)
  {
    fprintf(stderr, "Unexpected error: cannot cut the point in the inner target.\n");
    exit(EXIT_FAILURE);
  }

  // write_DLList(sol_neu_top_ddl, "SOL_NEU_TOP_DDL");
  // write_DLList(sol_neu_bottom_ddl, "SOL_NEU_BOTTOM_DDL");
  // write_DLList(sol_neu_left_ddl, "SOL_NEU_LEFT_DDL");
  // write_DLList(sol_neu_right_ddl, "SOL_NEU_RIGHT_DDL");

  Curve* sol_neu_top_curve=convert_ddl_to_curve(sol_neu_top_ddl);
  Curve* sol_neu_bottom_curve=convert_ddl_to_curve(sol_neu_bottom_ddl);
  Curve* sol_neu_left_curve=convert_ddl_to_curve(sol_neu_left_ddl);
  Curve* sol_neu_right_curve=convert_ddl_to_curve(sol_neu_right_ddl);

  //6.2 build the normal distribution
  int nbottom = sol_neu_bottom_curve->n_point;
  int nleft_sol=7;
  
  double* distrb_b=malloc((nbottom)*sizeof(double));
  double* distrb_l=malloc((nleft_sol)*sizeof(double));

  double len_tmp = total_length_curve(sol_neu_bottom_curve);

  for(int i=0; i<nbottom; i++)
  {
    distrb_b[i]=length_curve(sol_neu_bottom_curve, i+1)/len_tmp;
  }

  for(int i=0; i<nleft_sol; i++)
  {
    distrb_l[i]=i*(1.0/(nleft_sol-1));
  }

  distrb_b[0]=0.0;
  distrb_l[0]=0.0;
  distrb_b[nbottom-1]=1.0;
  distrb_l[nleft_sol-1]=1.0;
  
  write_array(distrb_b, nbottom, "distrb_b");
  write_array(distrb_l, nleft_sol, "distrb_l");
  //6.3 build the sol_neu 2d grid

  TwoDimGrid* sol_neu_2dgrid=create_2Dgrid_poloidal_major(nbottom, nleft_sol);

  generate_2Dgrid_default_TFI(sol_neu_2dgrid,
                              sol_neu_bottom_curve, sol_neu_top_curve, sol_neu_left_curve, sol_neu_right_curve,
                              distrb_b, nbottom,
                              distrb_b, nbottom,
                              distrb_l, nleft_sol,
                              distrb_l, nleft_sol);

  write_2Dgrid(sol_neu_2dgrid, "SOL_NEU_2DBASE");

  optimized_neu_2Dgrid(sol_neu_2dgrid);

  write_2Dgrid(sol_neu_2dgrid, "SOL_NEU_2DBASE_OPT");
  
  //6.4 free memmory
  free(distrb_b);
  free(distrb_l);

  free_DLList(sol_neu_top_ddl);
  free_DLList(sol_neu_bottom_ddl);
  free_DLList(sol_neu_left_ddl);
  free_DLList(sol_neu_right_ddl);
  
  free_curve(sol_neu_top_curve);
  free_curve(sol_neu_bottom_curve);
  free_curve(sol_neu_left_curve);
  free_curve(sol_neu_right_curve);

  free_2Dgrid(sol_neu_2dgrid);


/***********************************************
*   7. Expand the PFR neutral grid
***********************************************/
  //7.1 build the four boundaries
  DLListNode* pfr_neu_top_ddl = load_DLList_from_file("PFR_neu_bnd");
  
  r_tmp = get_x_2Dgrid(pfr2dgrid, 0, pfr2dgrid->nrad-1);
  z_tmp = get_y_2Dgrid(pfr2dgrid, 0, pfr2dgrid->nrad-1);
  DLListNode* pfr_neu_bottom_ddl = create_DLListNode(r_tmp, z_tmp);
  tail_tmp = pfr_neu_bottom_ddl;

  for(int i = 1; i<pfr2dgrid->npol;i++)
  {
    double r_tmp = get_x_2Dgrid(pfr2dgrid, i, pfr2dgrid->nrad-1);
    double z_tmp = get_y_2Dgrid(pfr2dgrid, i, pfr2dgrid->nrad-1);
    add_DLListnode_at_tail(&tail_tmp, r_tmp, z_tmp);
  }

  DLListNode* pfr_neu_left_ddl = load_DLList_from_file("inner_targetcurve");
  DLListNode* pfr_neu_right_ddl = load_DLList_from_file("outer_targetcurve");

  r_tmp = get_x_2Dgrid(pfr2dgrid,  pfr2dgrid->npol-1,pfr2dgrid->nrad-1);
  z_tmp = get_y_2Dgrid(pfr2dgrid,  pfr2dgrid->npol-1,pfr2dgrid->nrad-1);

  if(insert_point_for_DLList(pfr_neu_right_ddl, r_tmp, z_tmp))
  {
    fprintf(stderr, "Unexpected error: the point is not in the outer target.\n");
    exit(EXIT_FAILURE);
  }

  if(cut_DLList_after_point(pfr_neu_right_ddl, r_tmp, z_tmp)==0)
  {
    fprintf(stderr, "Unexpected error: cannot cut the point in the outer target.\n");
    exit(EXIT_FAILURE);
  }

  r_tmp = get_x_2Dgrid(pfr2dgrid, 0, pfr2dgrid->nrad-1);
  z_tmp = get_y_2Dgrid(pfr2dgrid, 0, pfr2dgrid->nrad-1);

  if(insert_point_for_DLList(pfr_neu_left_ddl, r_tmp, z_tmp))
  {
    fprintf(stderr, "Unexpected error: the point is not in the inner target.\n");
    exit(EXIT_FAILURE);
  }

  if(cut_DLList_after_point(pfr_neu_left_ddl, r_tmp, z_tmp)==0)
  {
    fprintf(stderr, "Unexpected error: cannot cut the point in the inner target.\n");
    exit(EXIT_FAILURE);
  }

  reverse_DLList(&pfr_neu_left_ddl);
  reverse_DLList(&pfr_neu_right_ddl);
  
  write_DLList(pfr_neu_top_ddl, "PFR_NEU_TOP_DDL");
  write_DLList(pfr_neu_bottom_ddl, "PFR_NEU_BOTTOM_DDL");
  write_DLList(pfr_neu_left_ddl, "PFR_NEU_LEFT_DDL");
  write_DLList(pfr_neu_right_ddl, "PFR_NEU_RIGHT_DDL");

  
  Curve* pfr_neu_top_curve=convert_ddl_to_curve(pfr_neu_top_ddl);
  Curve* pfr_neu_bottom_curve=convert_ddl_to_curve(pfr_neu_bottom_ddl);
  Curve* pfr_neu_left_curve=convert_ddl_to_curve(pfr_neu_left_ddl);
  Curve* pfr_neu_right_curve=convert_ddl_to_curve(pfr_neu_right_ddl);

  //7.2 build the normal distribution
  nbottom = pfr_neu_bottom_curve->n_point;
  nleft_sol=7;
  
  double* distrb_b_pfr=malloc((nbottom)*sizeof(double));
  double* distrb_l_pfr=malloc((nleft_sol)*sizeof(double));

  len_tmp = total_length_curve(pfr_neu_bottom_curve);

  for(int i=0; i<nbottom; i++)
  {
    distrb_b_pfr[i]=length_curve(pfr_neu_bottom_curve, i+1)/len_tmp;
  }

  for(int i=0; i<nleft_sol; i++)
  {
    distrb_l_pfr[i]=i*(1.0/(nleft_sol-1));
  }

  distrb_b_pfr[0]=0.0;
  distrb_l_pfr[0]=0.0;
  distrb_b_pfr[nbottom-1]=1.0;
  distrb_l_pfr[nleft_sol-1]=1.0;
  
  write_array(distrb_b, nbottom, "distrb_b_pfr");
  write_array(distrb_l, nleft_sol, "distrb_l_pfr");

  //7.3 build the sol_neu 2d grid

  TwoDimGrid* pfr_neu_2dgrid=create_2Dgrid_poloidal_major(nbottom, nleft_sol);

  generate_2Dgrid_default_TFI(pfr_neu_2dgrid,
                              pfr_neu_bottom_curve,pfr_neu_top_curve,pfr_neu_left_curve, pfr_neu_right_curve,
                              distrb_b_pfr, nbottom,
                              distrb_b_pfr, nbottom,
                              distrb_l_pfr, nleft_sol,
                              distrb_l_pfr, nleft_sol);

  write_2Dgrid(pfr_neu_2dgrid, "PFR_NEU_2DBASE");

  optimized_neu_2Dgrid(pfr_neu_2dgrid);

  write_2Dgrid(pfr_neu_2dgrid, "PFR_NEU_2DBASE_OPT");

  //7.4 free memmory
  free(distrb_b_pfr);
  free(distrb_l_pfr);

  free_DLList(pfr_neu_top_ddl);
  free_DLList(pfr_neu_bottom_ddl);
  free_DLList(pfr_neu_left_ddl);
  free_DLList(pfr_neu_right_ddl);
  
  free_curve(pfr_neu_top_curve);
  free_curve(pfr_neu_bottom_curve);
  free_curve(pfr_neu_left_curve);
  free_curve(pfr_neu_right_curve);

  free_2Dgrid(pfr_neu_2dgrid);

/**********************************************
*   8. Free space                             *
***********************************************/
  free(phi);
  // free_3Dgrid(sol3dgrid);
  // free_2Dgrid(grid_tmp1);
  // free_2Dgrid(grid_tmp2);
  // free_2Dgrid(grid_tmp3);

  // free_2Dgrid(sol2dgrid_exptgt);
  // free_2Dgrid(pfr2dgrid_exptgt);

  free_2Dgrid(sol2dgrid);
  free_2Dgrid(pfr2dgrid);

  // free_2Dgrid(core2dgrid);

  free_GridZone(solgz);
  free_GridZone(pfrgz);
  // free_GridZone(coregz);

  free_PolSegmsInfo(polseginfo);
  free_GridZoneInfo(&solgzinfo);
  free_GridZoneInfo(&pfrgzinfo);
  free_GridZoneInfo(&coregzinfo);

  free_dgtrg(trg);
  free_SepDistStr(sepdist);
  free_gradpsiline_default(gradpsilines);
  free_grad_psi(gradpsi);
  free_separatrix_default(sep);

  free_opoint(opoint);
  brk5_finalize(&brk45_data);
  free_interp1d_function(interp);
  free_2d_array(est_xpt);
  free_mag_field_torsys(&test_magfield);
  free_equilibrium(&dtt_example);

}

void EMC3_neu_3D_grid_generation_test()
{
  
/****************************
*  Read Necessary Inputs    *
*****************************/ 
  //inputs
  W3Config w3config;
  load_w3_config("w3.ini",&w3config);
  print_w3_config(&w3config);
 
  //Equilibirum
  Equilibrium dtt_example;
  init_equilibrium(&dtt_example);
  read_equilib_geqdsk(&dtt_example,w3config.file_config.geqdsk_file);
  correct_direction_lower_divertor(&dtt_example);
  print_equilibrium(&dtt_example);
  
  //Magnetic
  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);

  //DivGeo trg

  char* trgname=w3config.file_config.divgeo_trg_file;
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, trgname);

/************************************************
*  Build the 3D tracer                          *
************************************************/ 
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
  brk45_solver.initialize(&brk45_data);

/*****************
*  Read inputs   *
******************/ 

  //load SOL plasma 3D grid
  ThreeDimGrid* grid3d_plasma_SOL=load_EMC3_format_3Dgrid_from_file("grid3D_SOL.dat");
  ThreeDimGrid* grid3d_plasma_PFR=load_EMC3_format_3Dgrid_from_file("grid3D_PFR.dat");

  //load inner target
  DLListNode* sol_neu_left_ddl = load_DLList_from_file("inner_targetcurve");

  //load inner target
  DLListNode* sol_neu_right_ddl = load_DLList_from_file("outer_targetcurve");

  //load neutral boundary for SOL region
  DLListNode* sol_neu_top_ddl = load_DLList_from_file("SOL_neu_bnd");

  //load neutral boundary for PFR region
  DLListNode* pfr_neu_top_ddl = load_DLList_from_file("PFR_neu_bnd");

  //left bnd of neutral 3d grid. also the radial size
  //points number = cells number + 1
  int nleft_sol=w3config.grid2d_config.neu_expand_number[0]+1;
  double* distrb_l_sol=malloc((nleft_sol)*sizeof(double));

  for(int i=0; i<nleft_sol; i++)
  {
    distrb_l_sol[i]=i*(1.0/(nleft_sol-1));
  }

  distrb_l_sol[0]=0.0;
  distrb_l_sol[nleft_sol-1]=1.0;


  //points number = cells number + 1
  int nleft_pfr=w3config.grid2d_config.neu_expand_number[1]+1;
  double* distrb_l_pfr=malloc((nleft_pfr)*sizeof(double));

  for(int i=0; i<nleft_pfr; i++)
  {
    distrb_l_pfr[i]=i*(1.0/(nleft_pfr-1));
  }

  distrb_l_pfr[0]=0.0;
  distrb_l_pfr[nleft_pfr-1]=1.0;

  
  //For SOL neutral 3D GRID
  ThreeDimGrid* grid3d_neu_SOL = create_3Dgrid_radial_major(grid3d_plasma_SOL->npol, nleft_sol, grid3d_plasma_SOL->ntor);

  generate_EMC3_neutral_3Dgrid_TFI(grid3d_neu_SOL, grid3d_plasma_SOL,
                                  sol_neu_left_ddl, sol_neu_right_ddl, sol_neu_top_ddl,
                                  nleft_sol, distrb_l_sol, &ode_func, &brk45_solver);

  write_EMC3_3Dgrid_to_EMC3_format(grid3d_neu_SOL, "grid3D_NEU_SOL.dat",false,false);


  //For PFR neutral 3D GRID
  ThreeDimGrid* grid3d_neu_PFR = create_3Dgrid_radial_major(grid3d_plasma_PFR->npol, nleft_sol, grid3d_plasma_PFR->ntor);

  generate_EMC3_neutral_3Dgrid_TFI(grid3d_neu_PFR, grid3d_plasma_PFR,
                                   sol_neu_left_ddl, sol_neu_right_ddl, pfr_neu_top_ddl,
                                   nleft_pfr, distrb_l_pfr, &ode_func, &brk45_solver);
   
  write_EMC3_3Dgrid_to_EMC3_format(grid3d_neu_PFR, "grid3D_NEU_PFR.dat",false,false);


  /***************************************
  *   Merge plasma and neutral 3D grids  *
  ****************************************/ 
  ThreeDimGrid* grid3d_all_SOL=merge_3Dgrids_radial(grid3d_plasma_SOL,grid3d_neu_SOL);
  write_EMC3_3Dgrid_to_EMC3_format(grid3d_all_SOL, "grid3D_all_SOL.dat",false,false);

  ThreeDimGrid* grid3d_all_PFR=merge_3Dgrids_radial(grid3d_plasma_PFR,grid3d_neu_PFR);
  write_EMC3_3Dgrid_to_EMC3_format(grid3d_all_PFR, "grid3D_all_PFR.dat",false,false);


  // /***********************
  // *   Write BFIELD file  *
  // ************************/
  //   write_BFIELD_file_default(grid3d_all_SOL,&test_magfield,"BFIELD_SOL");
  //   write_BFIELD_file_default(grid3d_all_PFR,&test_magfield,"BFIELD_PFR");
  //   write_BFIELD_file_default(grid3d_plasma_core,&test_magfield,"BFIELD_CORE");

  // /**************************
  // *   Write PLATE_MAG file  *
  // ***************************/
  //   write_PLATE_MAG_file_test(grid3d_all_SOL, nleft_sol-1, 0, 1, "plate_SOL");
  //   write_PLATE_MAG_file_test(grid3d_all_PFR, nleft_sol-1, 0, 1, "plate_PFR");



  /****************************************
  *   Radial Mapping Check and Correction
  *****************************************/

  ThreeDimGrid* grid3d_plasma_core=load_EMC3_format_3Dgrid_from_file("grid3D_CORE.dat");
  
  close_core_3Dgrid_test(grid3d_plasma_core);

  int pol_inner=trg->nptseg[1]-1; //inner leg cell number, not point number
  int ntor=grid3d_plasma_core->ntor;

  int p1_s=0;
  int p1_e=grid3d_plasma_core->npol-1;
  int idx_r1=0;

  int p2_s=pol_inner+ntor/2;
  int p2_e=p2_s+p1_e-p1_s;
  int idx_r2=0;
  radial_mapping_check_test(grid3d_plasma_core,idx_r1,p1_s,p1_e,
                            grid3d_all_SOL,    idx_r2,p2_s,p2_e);
   
  p1_s=0;
  p1_e=pol_inner+ntor/2;
  p2_s=0;
  p2_e=pol_inner+ntor/2;
  idx_r1=0;
  idx_r1=0;
  radial_mapping_check_test(grid3d_all_SOL, idx_r1, p1_s, p1_e,
                            grid3d_all_PFR, idx_r2, p2_s, p2_e);

  p1_s=p1_e+grid3d_plasma_core->npol-1+1;
  p1_e=grid3d_all_SOL->npol-1;
    
  p2_s=p2_e+1;
  p2_e=grid3d_all_PFR->npol-1;
  radial_mapping_check_test(grid3d_all_SOL, idx_r1, p1_s, p1_e,
                            grid3d_all_PFR, idx_r2, p2_s, p2_e);

  /***************************************************
  *  Poloidal Expansion at inner and outer portions  *
  ****************************************************/
  int nsteps=10;
  poloidal_extend_for_toroidal_mapping_test(grid3d_all_SOL, nsteps, &ode_func, &brk45_solver);
  poloidal_extend_for_toroidal_mapping_test(grid3d_all_PFR, nsteps, &ode_func, &brk45_solver);


  write_EMC3_3Dgrid_to_EMC3_format(grid3d_all_SOL, "grid3D_all_SOL_polext.dat",false,false);
  write_EMC3_3Dgrid_to_EMC3_format(grid3d_all_PFR, "grid3D_all_PFR_polext.dat",false,false);


  /**********************************
  *   FOR T.LUNT EMC3 ORDER         *
  ***********************************/

  write_EMC3_3Dgrid_to_EMC3_format(grid3d_all_SOL, "grid3D_all_SOL_order.dat",true,false);
  write_EMC3_3Dgrid_to_EMC3_format(grid3d_all_PFR, "grid3D_all_PFR_order.dat",true,false);
  write_EMC3_3Dgrid_to_EMC3_format(grid3d_plasma_core, "grid3D_plasma_CORE_order.dat",true,true);


  /************************
  *    READ 3DGIRD        *
  *************************/

  ThreeDimGrid* grid3d_all_SOL_order=load_EMC3_format_3Dgrid_from_file("grid3D_all_SOL_order.dat");
  ThreeDimGrid* grid3d_all_PFR_order=load_EMC3_format_3Dgrid_from_file("grid3D_all_PFR_order.dat");
  ThreeDimGrid* grid3d_plasma_core_order=load_EMC3_format_3Dgrid_from_file("grid3D_plasma_CORE_order.dat");

 /**************************
  *    WRITE PLATE_MAG     *
  **************************/

  //OLD METHOD
  // write_PLATE_MAG_file_test(grid3d_all_SOL_order, nleft_sol-1, 1, 1, "plate_SOL_order");
  // write_PLATE_MAG_file_test(grid3d_all_PFR_order, nleft_sol-1, 1, 2, "plate_PFR_order");

  //NEW METHOD FOR SOL
  PlateCells* platecell_SOL=generate_platecells(grid3d_all_SOL_order, 1);
  update_platecells_targets(platecell_SOL, grid3d_all_SOL_order);
  update_platecells_radneu(platecell_SOL,w3config.grid2d_config.neu_expand_number[0]);

  DLListNode* limiter_head=load_DLList_from_file("antenna.s_e");
  double limiter_star_phi=0.0;
  double limiter_end_phi=24.0;
  update_platecells_limiter(platecell_SOL, grid3d_all_SOL_order,
                            limiter_head, limiter_star_phi, limiter_end_phi);
  

  //NEW METHOD FOR PFR
  PlateCells* platecell_PFR=generate_platecells(grid3d_all_PFR_order, 2);
  update_platecells_targets(platecell_PFR, grid3d_all_PFR_order);
  update_platecells_radneu(platecell_PFR,w3config.grid2d_config.neu_expand_number[1]);

  
  write_platecells(platecell_SOL,"plate_SOL_order_new");
  write_platecells(platecell_PFR,"plate_PFR_order_new");


  /*************************************************
  *  ONLY FOR TOROIDAL-AXIS-SYSMMETRY LIMITER      *
  **************************************************/
  //WRITE ANTENAA S_E SURFACE
  write_tor_axi_sym(grid3d_all_SOL_order, limiter_head, true, "ADD_antenna.s_e");




  free_DLList(limiter_head);
  free_platecells(platecell_SOL);

  /*************************************************
  *  ONLY FOR NON-TOROIDAL-AXIS-SYSMMETRY LIMITER  *
  **************************************************/

  //WRITE ANTENAA TOROIDAL SURFACE 
  // DLListNode* limiter_inner_bnd=load_DLList_from_file("antenna.inner.bnd.modify");
  // DLListNode* limiter_outer_bnd=load_DLList_from_file("antenna.outer.bnd.modify");

  // double phi_delta=1.0E-4;

  // write_non_tor_axi_sym(grid3d_all_SOL_order, limiter_head, false,
  //                       limiter_star_phi, limiter_end_phi, "ADD_antenna.s_e.modify");

  // write_non_tor_axi_sym_usr(limiter_head, false, limiter_star_phi,limiter_end_phi,
  //                           0.5,  "ADD_antenna.s_e.refine");

  // write_one_toroidal_surface(limiter_outer_bnd,limiter_star_phi-phi_delta,
  //                            limiter_inner_bnd,limiter_star_phi, 
  //                            false, "ADD_antenna.start.modify");
  
  // write_one_toroidal_surface(limiter_inner_bnd, limiter_end_phi,
  //                            limiter_outer_bnd, limiter_end_phi+phi_delta,
  //                            false, "ADD_antenna.end.modify");

  // free_DLList(limiter_inner_bnd);
  // free_DLList(limiter_outer_bnd);
  

 /***********************
  *    WRITE BFIELD     *
  ***********************/

  write_BFIELD_file_default(grid3d_all_SOL_order,&test_magfield,"BFIELD_SOL_order");
  write_BFIELD_file_default(grid3d_all_PFR_order,&test_magfield,"BFIELD_PFR_order");
  write_BFIELD_file_default(grid3d_plasma_core_order,&test_magfield,"BFIELD_CORE_order");


/***************************
*  Write OMP and Target RZ *
****************************/
  write_RZ_along_radial_test(grid3d_plasma_core_order, 43,0,"RZ_OMP_CORE");
  write_RZ_along_radial_test(grid3d_all_SOL_order, 99,0,"RZ_OMP_SOL");
  write_RZ_along_radial_test(grid3d_all_SOL_order, 12,0,"RZ_OT_SOL");
  write_RZ_along_radial_test(grid3d_all_SOL_order, 2,0,"RZ_IT_SOL");

/******************************
*  Write Additional surfaces  *
*******************************/
  int PFR_nr=grid3d_all_PFR_order->nrad;
  int SOL_nr=grid3d_all_SOL_order->nrad;


  int PFR_np=grid3d_all_PFR_order->npol;
  int SOL_np=grid3d_all_SOL_order->npol;
  int nt=grid3d_all_PFR_order->ntor;

  DLListNode* add_outer=create_DLListNode(get_r_3Dgrid(grid3d_all_PFR_order,nt,PFR_nr-1,0),
                                          get_z_3Dgrid(grid3d_all_PFR_order,nt,PFR_nr-1,0));

  DLListNode* add_inner=create_DLListNode(get_r_3Dgrid(grid3d_all_PFR_order,PFR_np-1-1,PFR_nr-1,0),
                                          get_z_3Dgrid(grid3d_all_PFR_order,PFR_np-1-1,PFR_nr-1,0));

  DLListNode* add_outer_tail=add_outer;
  DLListNode* add_inner_tail=add_inner;
  for(int i=1;i<PFR_nr;i++)
  {
    add_DLListnode_at_tail(&add_outer_tail,get_r_3Dgrid(grid3d_all_PFR_order,nt,PFR_nr-1-i,0),
                                           get_z_3Dgrid(grid3d_all_PFR_order,nt,PFR_nr-1-i,0));
    add_DLListnode_at_tail(&add_inner_tail,get_r_3Dgrid(grid3d_all_PFR_order,PFR_np-1-1,PFR_nr-1-i,0),
                                           get_z_3Dgrid(grid3d_all_PFR_order,PFR_np-1-1,PFR_nr-1-i,0));
  }

  for(int i=1;i<SOL_nr;i++)
  {
    add_DLListnode_at_tail(&add_outer_tail,get_r_3Dgrid(grid3d_all_SOL_order,nt,i,0),
                                           get_z_3Dgrid(grid3d_all_SOL_order,nt,i,0));
    add_DLListnode_at_tail(&add_inner_tail,get_r_3Dgrid(grid3d_all_SOL_order,SOL_np-1-1,i,0),
                                           get_z_3Dgrid(grid3d_all_SOL_order,SOL_np-1-1,i,0));
  }

  write_DLList(add_inner,"add_inner_ddl");
  write_DLList(add_outer,"add_outer_ddl");
  write_tor_axi_sym(grid3d_all_SOL_order, add_inner, true, "ADD_inner");
  write_tor_axi_sym(grid3d_all_SOL_order, add_outer, false, "ADD_outer");
  free_DLList(add_inner);
  free_DLList(add_outer);

/**************************
* Free parameters         *
**************************/ 
  free_dgtrg(trg);

  free_3Dgrid(grid3d_all_SOL_order);
  free_3Dgrid(grid3d_all_PFR_order);
  free_3Dgrid(grid3d_plasma_core_order);

  free(distrb_l_sol);
  free(distrb_l_pfr);

  free_3Dgrid(grid3d_all_SOL);
  free_3Dgrid(grid3d_all_PFR);
  free_3Dgrid(grid3d_plasma_core);

  free_3Dgrid(grid3d_plasma_SOL);
  free_3Dgrid(grid3d_neu_SOL);

  free_3Dgrid(grid3d_plasma_PFR);
  free_3Dgrid(grid3d_neu_PFR);

  free_DLList(sol_neu_left_ddl);
  free_DLList(sol_neu_right_ddl);
  free_DLList(sol_neu_top_ddl);
  free_DLList(pfr_neu_top_ddl);

  free_equilibrium(&dtt_example);
  free_mag_field_torsys(&test_magfield);
  brk5_finalize(&brk45_data);

}