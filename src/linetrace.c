#include "linetrace.h"

void initial_Bfield(Bfield_struct* Bfield)
{
  Bfield->nr = 0;
  Bfield->nz_RZ = 0;
  Bfield->nphi = 360/TOR_RESOLUTION + 1;
  Bfield->delta_phi = TOR_RESOLUTION;
  Bfield->r = NULL;
  Bfield->z_RZ = NULL;
  Bfield->Bfield_rzplane = NULL;
  
  Bfield->nx = 0;
  Bfield->ny = 0;
  Bfield->nz_XYZ = 0;
  Bfield->x = NULL;
  Bfield->y = NULL;
  Bfield->z_XYZ = NULL;
  Bfield->Bfield_xyz = NULL;
}

void create_Bfild(Bfield_struct* Bfield, const Equilibrium *equilib)
{
  Bfield->nr = equilib->nw;
  Bfield->nz_RZ = equilib->nh;
  Bfield->r = (double *)malloc(Bfield->nr * sizeof(double));

  for (int i=0; i<Bfield->nr; i++)
  {
    Bfield->r[i] = equilib->r[i];
  }

  Bfield->z_RZ = (double *)malloc(Bfield->nz_RZ * sizeof(double));

  for (int j=0; j<Bfield->nz_RZ;j++)
  {
    Bfield->z_RZ[j] = equilib->z[j];
  }
  
  Bfield->Bfield_rzplane = (double ***)malloc(Bfield->nr * sizeof(double **));
  Bfield->Bfield_rzplane[0] = (double **)malloc(Bfield->nr * Bfield->nz_RZ * sizeof(double *));
  Bfield->Bfield_rzplane[0][0] = (double *)malloc(Bfield->nr * Bfield->nz_RZ * 3 * sizeof(double));
  
  for (int i = 0; i < Bfield->nr; i++) 
  {
    Bfield->Bfield_rzplane[i] = Bfield->Bfield_rzplane[0] + i * Bfield->nz_RZ;
    for (int j = 0; j < Bfield->nz_RZ; j++) 
    {
      Bfield->Bfield_rzplane[i][j] = Bfield->Bfield_rzplane[0][0] + (i * Bfield->nz_RZ * 3) + (j * 3);
    }
  }
}

void psi_to_Bfield_rzplane(const Equilibrium *equilib, double ***Bfield_rzplane)
{
  derivation_2d(equilib->psi, equilib->r, equilib->nw, equilib->z, equilib->nh, Bfield_rzplane);
  // divide the field in R,Z direction by 2pi
  for(int i = 0; i<equilib->nw; i++)
  {
    for(int j=0; j<equilib->nh; j++)
    {
      Bfield_rzplane[i][j][0] = Bfield_rzplane[i][j][0] / PI / 2;
      Bfield_rzplane[i][j][1] = Bfield_rzplane[i][j][1] / PI / 2;
    }
  }
  double b0r0 = equilib->bcenter * equilib->rcenter;
  printf("debug: b0r0: %lf\n", b0r0);
  for(int i = 0; i<equilib->nw; i++)
  {
    for(int j=0; j<equilib->nh; j++)
    {
      if (equilib->r[i] < MIN_R) // To aviud unphysical values near R=0
      {
        Bfield_rzplane[i][j][2] = TOR_FIELD;
      }
      else
      {
        Bfield_rzplane[i][j][2] = b0r0 / equilib->r[i];
      }
      
    }
  }
}

void free_Bfield(Bfield_struct* Bfield)
{
  free(Bfield->Bfield_rzplane[0][0]);
  free(Bfield->Bfield_rzplane[0]);
  free(Bfield->Bfield_rzplane);
  free(Bfield->r);
  free(Bfield->z_RZ);
  free(Bfield->x);
  free(Bfield->y);
  free(Bfield->z_XYZ);
  free(Bfield->Bfield_xyz);
}

void write_Bfield_rzplane(Bfield_struct* Bfield)
{
  const char* rfilename;
  const char* zfilename; 
  const char* tfilename;
  const char* polfilename;
  rfilename = "Bfield_in_R";
  zfilename = "Bfield_in_Z";
  tfilename = "Bfield_in_Tor";
  polfilename = "Bfield_in_Pol";
  
  FILE* rfile = fopen(rfilename, "w");
  for (int i=0; i<Bfield->nr;i++)
  {
    for (int j = 0; j<Bfield->nz_RZ;j++)
    {
      fprintf(rfile, "%lf  %lf  %lf\n", \
      Bfield->r[i], Bfield->z_RZ[j], Bfield->Bfield_rzplane[i][j][0]);
    }
  }
  fclose(rfile);
  printf("write the field in R direction in %s\n", rfilename);

  FILE* zfile = fopen(zfilename, "w");
  for (int i=0; i<Bfield->nr;i++)
  {
    for (int j = 0; j<Bfield->nz_RZ;j++)
    {
      fprintf(zfile, "%lf  %lf  %lf\n", \
      Bfield->r[i], Bfield->z_RZ[j], Bfield->Bfield_rzplane[i][j][1]);
    }
  }
  fclose(zfile);
  printf("write the field in Z direction in %s\n", zfilename);
 
  FILE* tfile = fopen(tfilename, "w");
  for (int i=0; i<Bfield->nr;i++)
  {
    for (int j = 0; j<Bfield->nz_RZ;j++)
    {
      fprintf(tfile, "%lf  %lf  %lf\n", \
      Bfield->r[i], Bfield->z_RZ[j], Bfield->Bfield_rzplane[i][j][2]);
    }
  }
  fclose(tfile);
  printf("write the field in Toroidal direction in %s\n", tfilename);
 
  FILE* polfile = fopen(polfilename, "w");
  double tmp;
  for (int i=0; i<Bfield->nr;i++)
  {
    for (int j = 0; j<Bfield->nz_RZ;j++)
    {
      tmp = Bfield->Bfield_rzplane[i][j][0] * Bfield->Bfield_rzplane[i][j][0]
           +Bfield->Bfield_rzplane[i][j][1] * Bfield->Bfield_rzplane[i][j][1];
      tmp = sqrt(tmp);
      fprintf(polfile, "%lf  %lf  %lf\n", \
      Bfield->r[i], Bfield->z_RZ[j], tmp);
    }
  }
  fclose(polfile);
  printf("write the field in Poloidal direction in %s\n", polfilename);
}
