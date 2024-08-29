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

  find_Xpoint(&dtt_example, w3_input.xpt_estimation);
  

  
  
  // for(int j = 0; j < 5; j++)
  // {
  //   printf("z= %f\n", dtt_example.z[j]);
  // }
  
  //   for(int j = dtt_example.nh - 5; j <  dtt_example.nh; j++)
  // {
  //   printf("z= %f\n", dtt_example.z[j]);
  // }
 
  // test magnetic surface line 
  printf("psi is %lf\n",dtt_example.psi[120][120]);
  printf("psi is %lf\n",dtt_example.psi[121][120]);
  printf("psi is %lf\n",dtt_example.psi[120][121]);
  printf("psi is %lf\n",dtt_example.psi[120][121]);
  int test;
  test = calc_surface_line(&dtt_example,120,120,-0.13,dtt_example.nw, dtt_example.nh);
  printf("%d",test);
  free_equilibrium(&dtt_example);

  

  return 0;
}
