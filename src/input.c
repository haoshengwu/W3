#include <stdio.h>
#include <string.h>
#include "input.h"

// Here are the default values of the input parameters for building the mesh
// At the begining, we assume is single null. It the future, it should suport
// different configurations.


// Current use DTT values, will be updated in the future.
void init_inputpara(InputPara* input){

  strcpy(input -> equilibrium_file,"example.geqdsk");
  strcpy(input -> topology_type,"sinlge_null");
  //for DTT
  input->xpt_estimation[0] = 1.8617;  //  R coordinate
  input->xpt_estimation[1] = -1.1622; //  Z coordinate
  //for SPARC
  //input->xpt_estimation[0] = 1.5417;  //  R coordinate
  //input->xpt_estimation[1] = -1.1167; //  Z coordinate
}

void print_inputpara(InputPara* input){

  printf("Equilibirum file name: %s\n", input -> equilibrium_file);
  printf("The topology of equilibrium file is: %s\n", input -> topology_type);
  
}

