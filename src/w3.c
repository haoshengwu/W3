/* This is the main program of W3, which is a small project 
to crete the 3D emc3 mesh. The main aims are:
1. Create 1D plasma mesh, which can be used by ReMkit1D.
2. Create 2D plasma mesh, which can be used by SOLPS-ITER. 
3. Create 3D mesh, which can be used by EMC3-EIRENE.
Contact: haosheng.wu@polito.it
*/

#include <stdio.h>
#include <stdlib.h>

#include "test.h"

int main(){
  
  printf("Welcome to W3!\n");
  
  printf("If there is problem, please conatact haosheng.wu@polito.it\n");
  

//**********Test one line of the separatrix tracing**********************/
  // DLListNode* line_list=NULL;
  // cal_separatrix_line(&dtt_example, xp, 0,&(line_list));
  // write_DLList(line_list,"sep0");
  // free_DLList(line_list);
//**********************************************************************/


/*
test magnetic field line calculation
*/

//   MagFieldTorSys test_magfield;
//   init_mag_field_torsys(&test_magfield);
//   char* method = "central_4th";
//   calc_mag_field_torsys(&dtt_example, &test_magfield, method);
//   write_mag_field_torsys(&test_magfield);
//   printf("debug1\n");
//   printf("debug2\n");

// /*****************************************************
// test new structure for euler  tracing
// ******************************************************/
//   double direction[3]={-1.0,-1.0,-1.0};
  
//   ode_function ode_func = {
//     .ndim = 3,
//     .data = &test_magfield,
//     .rescale = direction,
//     .compute_f = ode_f_brz_torsys_bilinear,
//   };

//   ode_solver euler_solver =
//   {
//     .step_size = 1,
//     .solver_data = NULL,
//     .next_step = euler_next_step,
//     .initialize = euler_initialize,
//     .finalize = euler_finalize
//   };


//   int step = 360*10;
//   double *x1 = (double *)malloc((step+1) * sizeof(double));
//   x1[0] = 0.0;
//   double **line = allocate_2d_array(step+1,3);
//   line[0][0] = 2.75;
//   line[0][1] = 0.0;
//   line[0][2] = 0.0;

//   for(int i=1;i<step+1;i++)
//   {
//     x1[i] = x1[i-1] + euler_solver.step_size;
//     euler_solver.next_step(euler_solver.step_size, &(x1[i-1]), line[i-1], line[i], &test_magfield, &ode_func);
//   };
    
//   const char *filename1="euler_line_tracing";
//   FILE* file1 = fopen(filename1, "w");
//   for (int i=0; i<step+1; i++)
//   {
//       fprintf(file1, "%.12f %.12f %.12f\n", line[i][0],line[i][1],line[i][2]);
//   }
//   fclose(file1);
//   printf("write the tracing line in %s\n", filename1);

//   const char *filename2="euler_line_tracing_xyz";
//   FILE* file2 = fopen(filename2, "w");
//   for (int i=0; i<step+1; i++)
//   {
//       double x;
//       double y;
//       rphi_to_XY(line[i][0],line[i][2],&x,&y);
//       fprintf(file2, "%.12f %.12f %.12f\n", x,y,line[i][1]);
//   }
//   fclose(file2);
//   printf("write the tracing line in %s\n", filename2);

//   free(x1);
//   free_2d_array(line);

// /*****************************************************
// test new structure for brk5 tracing
// ******************************************************/
//   RKSolverData brk45_data;

//   ode_solver brk45_solver =
//   {
//     .step_size = 1,
//     .solver_data = &brk45_data,
//     .next_step = brk5_next_step,
//     .initialize = brk5_initialize,
//     .finalize = brk5_finalize
//   };

//   brk45_solver.initialize(&brk45_data);

//   int step2 =300;
//   double *x2 = (double *)malloc((step2+1) * sizeof(double));
//   x2[0] = 0.0;
//   double **line2 = allocate_2d_array(step2+1,3);
//   line2[0][0] = 2.75;
//   line2[0][1] = 0.0;
//   line2[0][2] = 0.0;

//   for(int i=1;i<step2+1;i++)
//   {
//     x2[i] = x2[i-1] + brk45_solver.step_size;
//     brk45_solver.next_step(brk45_solver.step_size, &(x2[i-1]), line2[i-1], line2[i], &brk45_data, &ode_func);
//   };
    
//   const char *filename3="brk5_line_tracing";
//   FILE* file3 = fopen(filename3, "w");
//   for (int i=0; i<step2+1; i++)
//   {
//       fprintf(file3, "%.12f %.12f %.12f\n", line2[i][0],line2[i][1],line2[i][2]);
//   }
//   fclose(file3);
//   printf("write the tracing line in %s\n", filename3);

//   const char *filename4="brk5_line_tracing_xyz";
//   FILE* file4 = fopen(filename4, "w");
//   for (int i=0; i<step2+1; i++)
//   {
//       double x;
//       double y;
//       rphi_to_XY(line2[i][0],line2[i][2],&x,&y);
//       fprintf(file4, "%.12f  %.12f  %.12f\n", x,y,line2[i][1]);
//   }
//   fclose(file4);
//   printf("write the tracing line in %s\n", filename4);

//   free(x2);
//   free_2d_array(line2);
//   brk45_solver.finalize(&brk45_data);



/*****************************************************
*  Verify the magnetic field line tracer
******************************************************/
  // tracer_test();
  // interpolator_test();
  //  line_tracer_test();

/*****************************************************
*  Verify the find_xpoint
******************************************************/
  // test_find_xpoint();

/*****************************************************
*  Verify the separatrix
******************************************************/
  // separatrix_test();
  // new_separatrix_test();
/*****************************************************
*  Verify Interpolation 1D
******************************************************/
  // interp1d_test();  

/*****************************************************
*  Verify read DG *.trg file
******************************************************/
  // read_trg_test();
  // get_target_curve_from_trg_test();

/*****************************************************
*  Verify target curve operations
******************************************************/
  // target_curve_test();

/*****************************************************
*  Verify grad psi
******************************************************/
  // grad_psi_test();


/*****************************************************
*  Verify divgeo
******************************************************/
  // divgeo_test();

/*****************************************************
*  Verify 2D carre grid generation 
******************************************************/
  meshgeneration_test();

  return 0;
}
