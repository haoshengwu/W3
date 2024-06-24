#include <stdio.h>
#include <stdlib.h>
#include "equilibrium.h"

void init_equilibrium(Equilibrium* equilib){

  equilib -> nw = 0;
  equilib -> nh = 0;
  equilib -> r = NULL;
  equilib -> z = NULL;
  equilib -> psi = NULL;

}

void print_equilibrium(const Equilibrium* equilib) {

  printf("Equilibrium size:\n");
  printf("nw: %i, * nh: %i\n", equilib -> nw, equilib -> nh);

}

//  Currently only store necessary paramters/
void read_equilib_geqdsk(Equilibrium* equilib, const char* geqdsk_file) {

  FILE* file = fopen(geqdsk_file,"r");

  if (file == NULL) {
    fprintf(stderr, "Error opening file %s\n", geqdsk_file);
    exit(1);  // Failed to read geqdsk equilibrium;
   }

  //  the length of title is 48 + 1 for /0;
  //  
  char title[49];
  int idum;
  double  rdim, zdim, rcentr, rleft, zmid;
  double  rmaxis, zmaxis, simag, sibry, bcentr;
  double  current,xdum;
  double  tmp;

  fscanf(file, "%48s %i %i %i", title, &idum, &equilib -> nw, &equilib -> nh);
  fscanf(file, "%lf %lf %lf %lf %lf", &rdim, &zdim, &rcentr, &rleft, &zmid);
  fscanf(file, "%lf %lf %lf %lf %lf", &rmaxis, &zmaxis, &simag, &sibry, &bcentr);
  fscanf(file, "%lf %lf %lf %lf %lf", &current, &simag, &xdum, &rmaxis, &xdum);
  fscanf(file, "%lf %lf %lf %lf %lf", &zmaxis, &xdum, &sibry, &xdum, &xdum);
 
  // skip fpol, pres, ffprim, pprime, total number is 4 * nw
  for(int i = 0; i < 4 * equilib -> nw; i++) {
    fscanf(file, "%lf", &tmp);
  }
 
  //  Allocate dynamic memmory and passing value for equilib -> psi
  equilib -> psi = (double**)malloc(equilib -> nw * sizeof(double*));

  for (int i = 0; i < equilib -> nw; i++) {

    equilib -> psi[i] = (double*)malloc(equilib -> nh * sizeof(double));

  }

  for (int j = 0; j < equilib -> nh; j++) {
    for (int i = 0; i < equilib -> nw; i++) {

      fscanf(file, "%lf", &equilib -> psi[i][j]);

    }
  }
 
  printf("finish reading geqdsk\n");
  // Allocate memmory and calculate r and z values.

  equilib -> r = (double*)malloc(equilib -> nw * sizeof(double));     
  equilib -> z = (double*)malloc(equilib -> nh * sizeof(double)); 

  double delta_r = rdim / (equilib -> nw - 1);
  for (int i = 0; i < equilib -> nw; i++) {

    equilib -> r[i] = i * delta_r + rleft;

  }

  double delta_z = zdim / (equilib -> nh - 1);
  double zleft = 2 * zmid - zdim;
  for (int j = 0; j < equilib -> nh; j++) {
  
    equilib -> z[j] = zleft + j * delta_z;
 
  }
  printf("r, z and psi is ok\n");

}


void free_equilibrium(Equilibrium* equilib) {
  for(int i = 0; i < equilib -> nw; i++) {
    free(equilib -> psi[i]);
    equilib -> psi[i] = NULL; 
  }

  free(equilib -> psi);
  equilib -> psi = NULL;
  
  free(equilib->r);
  free(equilib->z);
  equilib -> r = NULL;
  equilib -> z = NULL;
  
  equilib -> nw = 0;
  equilib -> nh = 0;
}
