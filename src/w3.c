/* This is the main program of W3, which is a small project 
to crete the 3D emc3 mesh. The main aims are:
1. Create 1D plasma mesh, which can be used by ReMkit1D.
2. Create 2D plasma mesh, which can be used by SOLPS-ITER. 
3. Create 3D mesh, which can be used by EMC3-EIRENE.
Contact: haosheng.wu@polito.it
*/

#include <stdio.h>
#include "input.h"

int main(){
  
  printf("Welcome to W3 project!\n");
  
  printf("If there is problem, please conatact haosheng.wu@polito.it\n");
  
  InputPara w3_input;
  
  init_inputpara(&w3_input);
  
  print_inputpara(&w3_input);

  return 0;
}
