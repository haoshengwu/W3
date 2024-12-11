/* This is the main program of W3, which is a small project 
to crete the 3D emc3 mesh. The main aims are:
1. Create 1D plasma mesh, which can be used by ReMkit1D.
2. Create 2D plasma mesh, which can be used by SOLPS-ITER. 
3. Create 3D mesh, which can be used by EMC3-EIRENE.
Contact: haosheng.wu@polito.it
*/

#include <stdio.h>
#include <stdlib.h>
#include "input.h"
#include "equilibrium.h"
#include "magneticsurface.h"
#include "datastructure.h"
#include "linetrace.h"
#include "magneticfield.h"
#include "ode.h"
#include "basemesh.h"

int main(){
  
  printf("Welcome to W3!\n");
  
  printf("If there is problem, please conatact haosheng.wu@polito.it\n");
  
  InputPara w3_input;
  
  init_inputpara(&w3_input);
  
  print_inputpara(&w3_input);

  Equilibrium dtt_example;
  
  init_equilibrium(&dtt_example);
 
  read_equilib_geqdsk(&dtt_example,w3_input.equilibrium_file);
  
  print_equilibrium(&dtt_example);
  XPointTest xp;
  xp = find_Xpoint(&dtt_example, w3_input.xpt_estimation);
  

  
  // for(int j = 0; j < 5; j++)
  // {
  //   printf("z= %f\n", dtt_example.z[j]);
  // }
  
  //   for(int j = dtt_example.nh - 5; j <  dtt_example.nh; j++)
  // {
  //   printf("z= %f\n", dtt_example.z[j]);
  // }
  // test module
  double value;
  double x = 1.534;
  double y = -1.125;
  value = get_psi_from_rz(&dtt_example,x, y);
  printf("psi at %lf %lf is %lf\n",x, y, value);

  // test magnetic surface line 
  //printf("psi is %lf\n",dtt_example.psi[149][126]);
  //printf("psi is %lf\n",dtt_example.psi[150][126]);
  //printf("psi is %lf\n",dtt_example.psi[149][127]);
  //printf("psi is %lf\n",dtt_example.psi[150][127]);
  
  //test separatrix
  // cal_separatrix_line(&dtt_example,xp,3);
  // int test;
  // test = calc_surface_line(&dtt_example,120,120,-0.135,dtt_example.nw, dtt_example.nh);
  // printf("%d",test);

  //test double linked list
  // double x1 = 1.0, y1 = 2.0;

  // DLListNode* ddl_ptr;
  // ddl_ptr = NULL;
  // ddl_ptr = create_DLListNode(x1,y1);
  // for (int i = 0; i < 10; i++)
  // {
  //   x1 = x1 + i;
  //   y1 = y1 + i;
  //   insert_DLList_at_head(&ddl_ptr,x1,y1);
  // }
  // print_DLList(ddl_ptr);
  

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

//**********Test one line of the separatrix tracing**********************/
  // DLListNode* line_list=NULL;
  // cal_separatrix_line(&dtt_example, xp, 0,&(line_list));
  // write_DDList(line_list,"sep0");
  // free_DLList(line_list);
//**********************************************************************/


/*
test magnetic field line calculation
*/

  MagFieldTorSys test_magfield;
  init_mag_field_torsys(&test_magfield);
  char* method = "central_2nd";
  calc_mag_field_torsys(&dtt_example, &test_magfield, method);
  write_mag_field_torsys(&test_magfield);
  printf("debug1\n");
  printf("debug2\n");


/*****************************************************
test new structure for euler  tracing
******************************************************/
  double direction[3]={-1.0,-1.0,-1.0};
  
  ode_function ode_func = {
    .ndim = 3,
    .data = &test_magfield,
    .rescale = direction,
    .compute_f = ode_f_brz_torsys_bilinear,
  };

  ode_solver euler_solver =
  {
    .step_size = 1,
    .solver_data = NULL,
    .next_step = euler_next_step,
    .initialize = euler_initialize,
    .finalize = euler_finalize
  };


  int step = 360*10;
  double *x1 = (double *)malloc((step+1) * sizeof(double));
  x1[0] = 0.0;
  double **line = allocate_2d_array(step+1,3);
  line[0][0] = 2.75;
  line[0][1] = 0.0;
  line[0][2] = 0.0;

  for(int i=1;i<step+1;i++)
  {
    x1[i] = x1[i-1] + euler_solver.step_size;
    euler_solver.next_step(euler_solver.step_size, &(x1[i-1]), line[i-1], line[i], &test_magfield, &ode_func);
  };
    
  const char *filename1="euler_line_tracing";
  FILE* file1 = fopen(filename1, "w");
  for (int i=0; i<step+1; i++)
  {
      fprintf(file1, "%.12f %.12f %.12f\n", line[i][0],line[i][1],line[i][2]);
  }
  fclose(file1);
  printf("write the tracing line in %s\n", filename1);

  const char *filename2="euler_line_tracing_xyz";
  FILE* file2 = fopen(filename2, "w");
  for (int i=0; i<step+1; i++)
  {
      double x;
      double y;
      rphi_to_XY(line[i][0],line[i][2],&x,&y);
      fprintf(file2, "%.12f %.12f %.12f\n", x,y,line[i][1]);
  }
  fclose(file2);
  printf("write the tracing line in %s\n", filename2);

  free(x1);
  free_2d_array(line);

/*****************************************************
test new structure for brk5 tracing
******************************************************/
  RKSolverData brk45_data;

  ode_solver brk45_solver =
  {
    .step_size = 1,
    .solver_data = &brk45_data,
    .next_step = brk5_next_step,
    .initialize = brk5_initialize,
    .finalize = brk5_finalize
  };

  brk45_solver.initialize(&brk45_data);

  int step2 =300;
  double *x2 = (double *)malloc((step2+1) * sizeof(double));
  x2[0] = 0.0;
  double **line2 = allocate_2d_array(step2+1,3);
  line2[0][0] = 2.75;
  line2[0][1] = 0.0;
  line2[0][2] = 0.0;

  for(int i=1;i<step2+1;i++)
  {
    x2[i] = x2[i-1] + brk45_solver.step_size;
    brk45_solver.next_step(brk45_solver.step_size, &(x2[i-1]), line2[i-1], line2[i], &brk45_data, &ode_func);
  };
    
  const char *filename3="brk5_line_tracing";
  FILE* file3 = fopen(filename3, "w");
  for (int i=0; i<step2+1; i++)
  {
      fprintf(file3, "%.12f %.12f %.12f\n", line2[i][0],line2[i][1],line2[i][2]);
  }
  fclose(file3);
  printf("write the tracing line in %s\n", filename3);

  const char *filename4="brk5_line_tracing_xyz";
  FILE* file4 = fopen(filename4, "w");
  for (int i=0; i<step2+1; i++)
  {
      double x;
      double y;
      rphi_to_XY(line2[i][0],line2[i][2],&x,&y);
      fprintf(file4, "%.12f  %.12f  %.12f\n", x,y,line2[i][1]);
  }
  fclose(file4);
  printf("write the tracing line in %s\n", filename4);

  free(x2);
  free_2d_array(line2);
  brk45_solver.finalize(&brk45_data);


/***************************************
 * Test the base mesh
 ***************************************/
//1. create the 2D trace solver
  double direction_2d[2] = {-1.0, -1.0};
  ode_function ode_func_2d = {
    .ndim = 2,
    .data = &test_magfield,
    .rescale = direction_2d,
    .compute_f = ode_f_brz_torsys_cubicherm,
  };

  RKSolverData brk45_data_2d;
  ode_solver brk45_solver_2d =
  {
    .step_size = 1,
    .solver_data = &brk45_data_2d,
    .next_step = brk5_next_step,
    .initialize = brk5_initialize,
    .finalize = brk5_finalize
  };
//2. Test 2D tracer: create 2 lines
  brk45_solver_2d.initialize(&brk45_data_2d);

  int step_prev= 1800;
  double *x_prev_2d = (double *)malloc((step_prev+1) * sizeof(double));
  x_prev_2d[0] = 0.0;
  double **line2d_prev = allocate_2d_array(step_prev+1,2);
  line2d_prev[0][0] = 2.05;
  line2d_prev[0][1] = -1.50;

//first line
  for(int i=1;i<step_prev+1;i++)
  {
    x_prev_2d[i] = x_prev_2d[i-1] + brk45_solver_2d.step_size;
    brk45_solver_2d.next_step(brk45_solver_2d.step_size, &(x_prev_2d[i-1]), line2d_prev[i-1], line2d_prev[i], &brk45_data_2d, &ode_func_2d);
  };
    
  const char *filename2d_prev="brk5_line_tracing_2d";
  FILE* file2d_prev = fopen(filename2d_prev, "w");
  for (int i=0; i<step_prev+1; i++)
  {
      fprintf(file2d_prev, "%.12f %.12f\n", line2d_prev[i][0],line2d_prev[i][1]);
  }
  fclose(file2d_prev);
  printf("write the tracing line in %s\n", filename2d_prev);
  
//second line
  int step_curve = 1700;
  double *x_2d = (double *)malloc((step_curve+1) * sizeof(double));
  x_2d[0] = 0.0;
  double **line2d = allocate_2d_array(step_curve+1,2);
  line2d[0][0] = 2.1;
  line2d[0][1] = -1.55;


  for(int i=1;i<step_curve+1;i++)
  {
    x_2d[i] = x_2d[i-1] + brk45_solver_2d.step_size;
    brk45_solver_2d.next_step(brk45_solver_2d.step_size, &(x_2d[i-1]), line2d[i-1], line2d[i], &brk45_data_2d, &ode_func_2d);
  };
    
  const char *filename2d="brk5_line_tracing_2d_2";
  FILE* file2d = fopen(filename2d, "w");
  for (int i=0; i<step_curve+1; i++)
  {
      fprintf(file2d, "%.12f %.12f\n", line2d[i][0],line2d[i][1]);
  }
  fclose(file2d);
  printf("write the tracing line in %s\n", filename2d);

//3. Create the mesh tube


    int n_point = 201;                 // 总点数
    double del = 1.0 / (n_point - 1);

    int m = 11;                        // 前后各11个点
    int mid_points = n_point - 2 * m;  // 中间点的个数

    double normal[n_point];            // 存储归一化长度坐标

    double del1 = 0.02;
    double del2 = 0.02;

    double front_step = del1 / (m - 1);    // 前11个点的步长
    double back_step = del2 / (m - 1);     // 后11个点的步长
    double middle_step = (1 - del1 - del2) / (mid_points+1); // 中间区域步长

    // 前11个点：累积长度从 0.0 到 0.05
    normal[0] = 0.0;
    for (int i = 1; i < m; i++) {
        normal[i] = normal[i - 1] + front_step;
    }

    // 中间79个点：从前区域结束位置开始累加
    for (int i = m; i < n_point - m; i++) {
        normal[i] = normal[i - 1] + middle_step;
    }

    // 后11个点：从中间区域结束位置开始累加
    for (int i = n_point - m; i < n_point; i++) {
        normal[i] = normal[i - 1] + back_step;
    }

    // 修正最后一个点，确保精确为 1.0
    normal[n_point - 1] = 1.0;

    // 输出结果
    for (int i = 0; i < n_point; i++) {
        printf("normal[%d] = %.10f\n", i, normal[i]);
    }

  int n_prev_curve = step_prev+1;
  int n_curve = step_curve+1;

  printf("debug n_prev_curve: %d\n",n_prev_curve);
  printf("debug n_curve: %d\n",n_curve);
  
  double tot_length_prve = long_CARRE(line2d_prev,n_prev_curve);
  printf("total prev curve length: %f\n", tot_length_prve);

  double tot_length = long_CARRE(line2d,n_curve);
  printf("total curve length: %f\n", tot_length);

  double length_prev_points[n_point];
  double **prev_point_coord = allocate_2d_array(n_point,2);


   for(int i=0; i<n_point; i++)
   {
     printf("del: %f, i: %d, tot_length_prve: %.10f\n", del, i, tot_length_prve);
     length_prev_points[i] = normal[i] * tot_length_prve;
     printf("length_prev_points: %.10f\n",length_prev_points[i]);
     coord_CARRE(line2d_prev, n_prev_curve, length_prev_points[i], prev_point_coord[i]);
     printf("prev_point_coord: x: %.10f, y: %.10f\n",prev_point_coord[i][0],prev_point_coord[i][1]);
   }



  const char *prev_point_name="prev_point";
  FILE* prev_point = fopen(prev_point_name, "w");
  for (int i=0; i<n_point; i++)
  {
      fprintf(prev_point, "%.12f %.12f\n", prev_point_coord[i][0],prev_point_coord[i][1]);
  }
  fclose(prev_point);
  printf("write the tracing line in %s\n", prev_point_name);

  double guard_top = 0.05;
  double guard_end = 0.2;
  double pasmin = 1.0e-3;

  double length_points[n_point];
  // length_points[0] = 0;
  // length_points[n_point-1] = tot_length;

  double **point_coord = allocate_2d_array(n_point,2);
  //point_coord[0][0] = line2d_2[0][0];
  //point_coord[0][1] = line2d_2[0][1];
  //point_coord[n_point-1][0] = line2d_2[n_curve-1][0];
  //point_coord[n_point-1][1] = line2d_2[n_curve-1][1];
 
  CarreMeshTube tube1 =
  {
    .guard_top = guard_top,
    .guard_end = guard_end,
    .pasmin = pasmin,

    .n_prev_curve = n_prev_curve,
    .prev_curve = line2d_prev,

    .n_curve = n_curve,
    .curve = line2d,

    .n_point = n_point,
    .length_prev_points = length_prev_points,
    .prev_point_coord = prev_point_coord,

    .length_points = length_points,
    .point_coord = point_coord
  };


//4. Calculate the curve
  calc_points_CARRE(&tube1);

//5. output and free

  const char *filenametube="flux_tube";
  FILE* filetube = fopen(filenametube, "w");
  for (int i=0; i<n_point; i++)
  {
      fprintf(filetube, "%.12f %.12f\n", tube1.prev_point_coord[i][0],tube1.prev_point_coord[i][1]);
      fprintf(filetube, "%.12f %.12f\n", tube1.point_coord[i][0],tube1.point_coord[i][1]);
  }
  // for (int i=0; i<n_point; i++)
  // {
  //     fprintf(filetube, "%.12f %.12f\n", tube1.point_coord[i][0],tube1.point_coord[i][1]);
  // }
  fclose(filetube);
  printf("write the tracing line in %s\n", filenametube);

  brk45_solver_2d.finalize(&brk45_data_2d);
  free(x_2d);
  free(x_prev_2d);
  free_2d_array(line2d);
  free_2d_array(line2d_prev);
  free_2d_array(point_coord);
  free_2d_array(prev_point_coord);
  free_mag_field_torsys(&test_magfield);


  // Bfield_struct test_bfield;
  // initial_Bfield(&test_bfield);
  // create_Bfild(&test_bfield, &dtt_example);
  // psi_to_Bfield_rzplane(&dtt_example, test_bfield.Bfield_rzplane);
  // write_Bfield_rzplane(&test_bfield);




/*

  double delta_phi = 1;

  
  sparc value:
  int step = 1008;
  double r0 = 2.41;
  */
  /*
  ddt core
  int step = 504;
  double r0 = 2.8;
  */
  /*
  ddt SOL
  int step = 720;
  double r0 = 2.871;
  */
//   int step = 360;
//   double r0 = 2.871;
//   double phi0 = 0;
//   double z0 = 0.0;
//   double** line = allocate_2d_array(step+1,3);
  

//   euler_method(r0,z0,phi0, delta_phi, step, &test_bfield,1.0,1.0,line);

// // opposite direction
//   int step2 = 720;
//   double** line2 = allocate_2d_array(step2+1,3);
//   euler_method(r0,z0,phi0, -delta_phi, step2, &test_bfield,-1.0,-1.0,line2);

//   const char *filename3="euler_line_tracing_xyz_line2";
//   FILE* file3 = fopen(filename3, "w");
//   for (int i=0; i<step2+1; i++)
//   {
//       double x;
//       double y;
//       rphi_to_XY(line2[i][0],line2[i][1],&x,&y);
//       fprintf(file3, "%lf  %lf  %lf\n", x,y,line2[i][2]);
//   }
//   fclose(file3);
//   printf("write the tracing line in %s\n", filename3);


//   const char *filename1="euler_line_tracing";
//   FILE* file1 = fopen(filename1, "w");
//   for (int i=0; i<step+1; i++)
//   {
//       fprintf(file1, "%lf  %lf  %lf\n", line[i][0],line[i][1],line[i][2]);
//   }
//   fclose(file1);
//   printf("write the tracing line in %s\n", filename1);

//   const char *filename2="euler_line_tracing_xyz";
//   FILE* file2 = fopen(filename2, "w");
//   for (int i=0; i<step+1; i++)
//   {
//       double x;
//       double y;
//       rphi_to_XY(line[i][0],line[i][1],&x,&y);
//       fprintf(file2, "%lf  %lf  %lf\n", x,y,line[i][2]);
//   }
//   fclose(file2);
//   printf("write the tracing line in %s\n", filename2);

//   free_Bfield(&test_bfield);
//   free_2d_array(line);
//   free_2d_array(line2);
  
  free(xp);
  free_equilibrium(&dtt_example);  


  return 0;
}
