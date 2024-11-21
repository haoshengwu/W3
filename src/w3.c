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
  Bfield_struct test_bfield;
  initial_Bfield(&test_bfield);
  create_Bfild(&test_bfield, &dtt_example);
  psi_to_Bfield_rzplane(&dtt_example, test_bfield.Bfield_rzplane);
  write_Bfield_rzplane(&test_bfield);


/*
test magnetic field line calculation
*/
  double delta_phi = 0.5;

  /*
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
  int step = 720*3;
  double r0 = 1.490;
  double phi0 = 0;
  double z0 = 0.0;
  double** line = allocate_2d_array(step+1,3);
  for (int i = 0; i <= step; i++) 
  {
    line[i] = malloc(3 * sizeof(double));  // Allocate each row
  }
  euler_method(r0,z0,phi0, delta_phi, step, &test_bfield,1,1,line);

  const char *filename1="euler_line_tracing";
  FILE* file1 = fopen(filename1, "w");
  for (int i=0; i<step+1; i++)
  {
      fprintf(file1, "%lf  %lf  %lf\n", line[i][0],line[i][1],line[i][2]);
  }
  fclose(file1);
  printf("write the tracing line in %s\n", filename1);

  const char *filename2="euler_line_tracing_xyz";
  FILE* file2 = fopen(filename2, "w");
  for (int i=0; i<step+1; i++)
  {
      double x;
      double y;
      rphi_to_XY(line[i][0],line[i][1],&x,&y);
      fprintf(file2, "%lf  %lf  %lf\n", x,y,line[i][2]);
  }
  fclose(file2);
  printf("write the tracing line in %s\n", filename2);

  free_Bfield(&test_bfield);
  free_2d_array(line);
  free(xp);
  free_equilibrium(&dtt_example);  


  return 0;
}
