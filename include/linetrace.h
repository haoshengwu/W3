#ifndef LINETRACE_H
#define LINETRACE_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "equilibrium.h"
#include "calc.h"
#include "datastructure.h"

#define TOR_TURN 10
#define TOR_RESOLUTION 0.1
#define TOR_FIELD 200
#define MIN_R 0.1


typedef struct
{
  double rcenter;
  double bcenter;
  int nr;   // number or R point
  int nz_RZ;   // number of Z point
  int nphi; ;  // number of phi point
 
  double delta_phi; 
  double *r;
  double *z_RZ;
  double ***Bfield_rzplane;
  
  int nx;
  int ny;
  int nz_XYZ;
//  double ***x;
//  double ***y;
//  double ***z_XYZ;
//  double ****Bfield_xyz;
} Bfield_struct;

void initial_Bfield(Bfield_struct* Bfiled);

void create_Bfild(Bfield_struct* Bfield, const Equilibrium *equilib);

void psi_to_Bfield_rzplane(const Equilibrium *equilib, double ***Bfield_rzplane);


void free_Bfield(Bfield_struct* Bfiled);

void write_Bfield_rzplane(Bfield_struct* Bfield);

void euler_method(double r0, double z0, double phi0, double dphi, int step, 
                  Bfield_struct* Bfield, double pol_dir, double tor_dir,
                  double** output);
#endif