#ifndef __EQUILIBRIUM_H_
#define __EQUILIBRIUM_H_

//  Define the structure of equlibrium. Currently, it is assuemd the format is geqsk file.
//  
//  Related functions are also clear here.

#include <stdlib.h>
#include <string.h>

typedef struct {

  int  nw;    //  number of points along the R coordinate. 
  int  nh;    //  number of points along the Z coordinate.
  double *r, *z;    //  R and Z coordinates, dynamic arrays.
  double **psi;      //  Psi value, dynamic arrays.

}  Equilibrium;

void init_equilibrium(Equilibrium* equilib); //  initialize 

// read from geqsk file     
void read_equilib_geqdsk(Equilibrium* equilib, const char* geqdsk_file);

//  void write_equilib(Equilibrium* equilib, const char* filename) 
//  need to think the write format.

void print_equilibrium(const Equilibrium* equilib);  // print the size of equilibrium
void free_equilibrium(Equilibrium* equilib);         // free dynamic memory

#endif
