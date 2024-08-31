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
  double x = 1872.5/1000;
  double y = -1097.5/1000;
  value = get_psi_from_rz(&dtt_example,x, y);
  printf("psi at %lf %lf is %lf\n",x, y, value);

  // test magnetic surface line 
  printf("psi is %lf\n",dtt_example.psi[149][126]);
  printf("psi is %lf\n",dtt_example.psi[150][126]);
  printf("psi is %lf\n",dtt_example.psi[149][127]);
  printf("psi is %lf\n",dtt_example.psi[150][127]);
  
  //test separatrix
  // cal_separatrix_line(&dtt_example,xp,3);
  free(xp);
  // int test;
  // test = calc_surface_line(&dtt_example,120,120,-0.135,dtt_example.nw, dtt_example.nh);
  // printf("%d",test);
  free_equilibrium(&dtt_example);  

  //test double linked list
  double x1 = 1.0, y1 = 2.0;

  DLListNode* ddl_ptr;
  ddl_ptr = NULL;
  ddl_ptr = create_DLListNode(x1,y1);
  for (int i = 0; i < 10; i++)
  {
    x1 = x1 + i;
    y1 = y1 + i;
    ddl_ptr = insert_DLList_at_head(ddl_ptr,x1,y1);
  }
  print_DLList(ddl_ptr);
  free_DLList(ddl_ptr);

  return 0;
}
