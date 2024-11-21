#include "linetrace.h"

void initial_Bfield(Bfield_struct* Bfield)
{
  Bfield->rcenter = 0;
  Bfield->bcenter = 0;

  Bfield->nr = 0;
  Bfield->nz_RZ = 0;
  Bfield->nphi = 360/TOR_RESOLUTION + 1;
  Bfield->delta_phi = TOR_RESOLUTION;
  Bfield->r = NULL;
  Bfield->z_RZ = NULL;
  //Bfield_rzplane[i][j][0]:Br
  //Bfield_rzplane[i][j][1]:Bz
  //Bfield_rzplane[i][j][2]:Bphi
  Bfield->Bfield_rzplane = NULL;

  
  
  Bfield->nx = 0;
  Bfield->ny = 0;
  Bfield->nz_XYZ = 0;
//  Bfield->x = NULL;
//  Bfield->y = NULL;
//  Bfield->z_XYZ = NULL;
//  Bfield->Bfield_xyz = NULL;
}

void create_Bfild(Bfield_struct* Bfield, const Equilibrium *equilib)
{
  Bfield->bcenter = equilib->bcenter;
  Bfield->rcenter = equilib->rcenter;
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
  
  Bfield->Bfield_rzplane = allocate_3d_array(Bfield->nr,Bfield->nz_RZ,3);
  for(int i = 0; i<equilib->nw; i++)
  {
    for(int j=0; j<equilib->nh; j++)
    {
      for(int k=0; k<3; k++)
      {
        Bfield->Bfield_rzplane[i][j][k] = 10.0;
      }
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
      // due to the magnetic filed equation
      swap(Bfield_rzplane[i][j][0],Bfield_rzplane[i][j][1]);
      // Do not forget -1 for Br!!!
      if (equilib->r[i] < MIN_R) equilib->r[i] = MIN_R;
      Bfield_rzplane[i][j][0] = -Bfield_rzplane[i][j][0] / 2 /equilib->r[i];
      Bfield_rzplane[i][j][1] = Bfield_rzplane[i][j][1] / 2 /equilib->r[i];
    }
  }
  double b0r0 = equilib->bcenter * equilib->rcenter;
  printf("debug: b0r0: %lf\n", b0r0);
  for(int i = 0; i<equilib->nw; i++)
  {
    for(int j = 0; j<equilib->nh; j++)
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
//  free(Bfield->x);
//  free(Bfield->y);
//  free(Bfield->z_XYZ);
//  free(Bfield->Bfield_xyz);
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
  
  printf("***************************\n");
  printf("debug\n");
  printf("***************************\n");

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

  printf("***************************\n");
  printf("debug\n");
  printf("***************************\n");

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

  printf("***************************\n");
  printf("debug\n");
  printf("***************************\n");

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

/*
euler method to trace from the start point (r0,phi0,z0) along the magnetic filed in R Phi Z coordinate
start point: r0, phi0, z0
magnetic file info: Bfield
delta size: dphi
steps: step
rescale poloidal direction: pol_dir
rescale toroidal direction: tor_dir
results: output[0:step-1][0:2] 
output[:][0]: R
output[:][1]: phi
output[:][2]: Z
*/ 
void euler_method(double r0, double z0, double phi0, double dphi, int step, 
                  Bfield_struct* Bfield, double pol_dir, double tor_dir,
                  double** output)
{
  printf("Use euler method to trace line\n");
/*
!!! To do check the size of output satisfy the size: step
*/


//first point is the start point
  output[0][0] = r0;
  output[0][1] = phi0;
  output[0][2] = z0;

//get temporary values for Br(0:nh-1,0:nw-1) and Bz (0:nh-1,0:nw-1)
  double **br_tmp = allocate_2d_array(Bfield->nr,Bfield->nz_RZ);
  double **bz_tmp = allocate_2d_array(Bfield->nr,Bfield->nz_RZ);
  for(int i=0;i<Bfield->nr;i++)
  {
    for(int j=0;j<Bfield->nz_RZ;j++)
    {
      br_tmp[i][j] = Bfield->Bfield_rzplane[i][j][0];
      bz_tmp[i][j] = Bfield->Bfield_rzplane[i][j][1];
    }
  }
  // coeffcients in each cell, so the number of cell is nr-1 and nz-1, not same with point number.
  double ***br_bilinear_coeff = allocate_3d_array(Bfield->nr-1,Bfield->nz_RZ-1,4);
  double ***bz_bilinear_coeff = allocate_3d_array(Bfield->nr-1,Bfield->nz_RZ-1,4);

  bilinear_coeff(br_tmp, Bfield->r, Bfield->z_RZ, Bfield->nr, Bfield->nz_RZ,br_bilinear_coeff);
  bilinear_coeff(bz_tmp, Bfield->r, Bfield->z_RZ, Bfield->nr, Bfield->nz_RZ,bz_bilinear_coeff);

  double b0_tmp = Bfield->bcenter;
  double r0_tmp = Bfield->rcenter;


//
  for (int i=1; i<step+1; i++)
  {
    printf("The ith point: %d\n", i);
    //use previous piont as a temporary value
    double r_tmp = output[i-1][0];
    double z_tmp = output[i-1][2];
    printf("debug r_tmp: %lf, z_tmp: %lf\n",r_tmp,z_tmp);
    //ensure the cell number ii,jj
    int ii = ifind(r_tmp,Bfield->r,Bfield->nr,1);
    
    int jj = ifind(z_tmp,Bfield->z_RZ,Bfield->nz_RZ,1);
    printf("debug ii: %d, jj: %d\n",ii,jj);
    printf("debug br_tmp: %lf, bz_tmp: %lf\n",br_tmp[ii][jj],bz_tmp[ii][jj]);

    double xx = r_tmp;
    double zz = z_tmp;
    printf("debug xx: %lf, yy: %lf\n",xx,zz);

    double br_iijj = br_bilinear_coeff[ii][jj][0] + br_bilinear_coeff[ii][jj][1] * xx
                     + br_bilinear_coeff[ii][jj][2]*zz + br_bilinear_coeff[ii][jj][3] * xx * zz;
    br_iijj = br_iijj * pol_dir;
    
    double bz_iijj = bz_bilinear_coeff[ii][jj][0] + bz_bilinear_coeff[ii][jj][1] * xx
                     + bz_bilinear_coeff[ii][jj][2]*zz + bz_bilinear_coeff[ii][jj][3] * xx * zz;
    bz_iijj = bz_iijj * pol_dir;
    
    printf("debug br_iijj: %lf, bz_iijj: %lf\n",br_iijj,bz_iijj);

    double bphi_iijj = b0_tmp * r0_tmp / r_tmp;
    bphi_iijj = bphi_iijj * tor_dir;
    printf("debug bphi_iijj: %lf\n",bphi_iijj);
    
  // please do not forget use arc instead of phi
    double drad_tmp = dphi * PI / 180.0;
    
    output[i][0] = output[i-1][0] + r_tmp*drad_tmp*br_iijj/bphi_iijj*PI;
    output[i][1] = output[i-1][1] + dphi;
    output[i][2] = output[i-1][2] + r_tmp*drad_tmp*bz_iijj/bphi_iijj*PI;
  }

/*
!!! To do check whethe it is out of domain
*/
  free_2d_array(br_tmp);
  free_2d_array(bz_tmp);
  free_3d_array(br_bilinear_coeff);
  free_3d_array(bz_bilinear_coeff);
}