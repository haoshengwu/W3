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
  n_cut = cut_DLList_from_intersections(sep->line_list[num_itrsct], itrsct_r,itrsct_z);
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


  char* trgname="example.trg";
  DivGeoTrg* trg=create_dgtrg();
  int status=load_dgtrg_from_file(trg, trgname);

/***********************************************
*   Update trg regions(psi valuse)
***********************************************/
  for(int i=0; i<3; i++)
  {
    trg->regions[i]->level[0]=xpt_array[1].level;
  }

  write_sn_gridzoneinfo_from_dgtrg(trg, &dtt_example, sep, gradpsilines);
  


/***********************************************
*   Read input of GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_SOL");
  GridZoneInfo* pfrgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_PFR");
  GridZoneInfo* coregzinfo=load_GridZoneInfo_from_input("gridzoneinfo_CORE");

  print_GridZoneInfo(solgzinfo);
  print_GridZoneInfo(pfrgzinfo);
  print_GridZoneInfo(coregzinfo);


  write_polsegms_from_dgtrg(trg,"polseginfo1");
  PolSegmsInfo* polseginfo=read_PolSegmsInfo_from_file("polseginfo1");
  write_PolSegmsInfo(polseginfo,"polseginfo2");

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
*   Read input of GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_SOL");
  GridZoneInfo* pfrgzinfo=load_GridZoneInfo_from_input("gridzoneinfo_PFR");
  GridZoneInfo* coregzinfo=load_GridZoneInfo_from_input("gridzoneinfo_CORE");

/***********************************************
*   Read input of polsegminfo
***********************************************/
  PolSegmsInfo* polseginfo=read_PolSegmsInfo_from_file("polseginfo2");

  
/***********************************************
*   Create separatrix distribution
***********************************************/
  SepDistStr* sepdist=create_SepDistStr_from_sep(sep);
  update_sn_SepDistStr_from_GridZoneInfo(sepdist,solgzinfo);
  update_sn_SepDistStr_from_PolSegmsInfo(sepdist,polseginfo);
  update_SepDistStr_gridpoint_curve(sepdist);
  // int idx=sepdist->index[0];
  // write_curve("gridpoint_curve0", sepdist->edges[idx]->gridpoint_curve);
  // idx=sepdist->index[1];
  // write_curve("gridpoint_curve1", sepdist->edges[idx]->gridpoint_curve);
  // idx=sepdist->index[2];
  // write_curve("gridpoint_curve2", sepdist->edges[idx]->gridpoint_curve);

  // write_DLList(sepdist->edges[sepdist->index[0]]->head,"sepdist_line0");
  // write_DLList(sepdist->edges[sepdist->index[1]]->head,"sepdist_line1");
  // write_DLList(sepdist->edges[sepdist->index[2]]->head,"sepdist_line2");
  // write_DLList(sepdist->edges[sepdist->index[3]]->head,"sepdist_line3");


  GridZone* solgz=create_sn_GridZone(solgzinfo, sepdist);
  GridZone* pfrgz=create_sn_GridZone(pfrgzinfo, sepdist);
  GridZone* coregz=create_sn_GridZone(coregzinfo, sepdist);

  write_curve("sol_gz_c",solgz->first_bnd_curve);
  write_curve("sol_gz_gpc",solgz->first_gridpoint_curve);
  write_curve("pfr_gz_c",pfrgz->first_bnd_curve);
  write_curve("pfr_gz_gpc",pfrgz->first_gridpoint_curve);
  write_curve("core_gz_c",coregz->first_bnd_curve);
  write_curve("core_gz_gpc",coregz->first_gridpoint_curve);

  
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