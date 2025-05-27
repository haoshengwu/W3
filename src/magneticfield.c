#include "magneticfield.h"

#define NUM_DIFF_METHODS (sizeof(diff_2d_methods) / sizeof(diff_2d_methods[0]))
#define NUM_INTPL_METHODS (sizeof(intpl_2d_methods) / sizeof(intpl_2d_methods[0]))

DiffMethodEntry diff_2d_methods[] = {
  {"central_2nd", central_2nd_2d_diff},
  {"central_4th", central_4th_2d_diff},
// to do:
//  {"central_6th", central_6th_2d_diff},
//  {"central_8th", central_8th_2d_diff}
};

diff_2d_fun get_diff_method(const char *name)
{
  for (int i=0; i < NUM_DIFF_METHODS; i++)
  {
    if(strcmp(diff_2d_methods[i].name, name) == 0)
    {
      return diff_2d_methods[i].func;
    }
  }
  printf("The method: '%s' is not support\n", name);
  printf("ERROR: Please check the name of differential method\n");
  exit(1);
}

InterpolateMethodEntry intpl_2d_methods[] = {
  {"bilnear", bilenar2d2f},
  {"bicubic", bicubic2d2f},
  {"cubicherm", cubicherm2d2f},
//  {"cubichermite", cubherm_2d},
// to do:
};

intpl_2d_fun get_intpl_2d_method(const char *name)
{
  for (int i=0; i < NUM_INTPL_METHODS; i++)
  {
    if(strcmp(intpl_2d_methods[i].name, name) == 0)
    {
      return intpl_2d_methods[i].func;
    }
  }
  printf("The method: '%s' is not support\n", name);
  printf("ERROR: Please check the name of interpolated method\n");
  exit(1);
}

void init_mag_field_torsys(MagFieldTorSys *mag_field)
{
  mag_field->nr = 0;
  mag_field->nphi = 0;
  mag_field->nz = 0;
  mag_field->b0r0 = 0.0;
  mag_field->dphi = 0.0;
  mag_field->r = NULL;
  mag_field->z = NULL;
  mag_field->phi = NULL;
  mag_field->Brz = NULL;
  mag_field->dBrzdx = NULL;
  mag_field->dBrzdy = NULL;
  mag_field->d2Brzdxdy = NULL;
}

void free_mag_field_torsys(MagFieldTorSys *mag_field)
{
  free(mag_field->r);
  free(mag_field->z);
  free(mag_field->phi);
  free_3d_array(mag_field->Brz);
  free_3d_array(mag_field->dBrzdx);
  free_3d_array(mag_field->dBrzdy);
  free_3d_array(mag_field->d2Brzdxdy);
}

void calc_mag_field_torsys(Equilibrium *equ, MagFieldTorSys *mag_field, const char *method)
{
  diff_2d_fun selected_fun= get_diff_method(method);
  if (selected_fun == NULL)
  {
    printf("ERROE: Method: '%s' not found!\n", method);
    exit(1);
  }
  mag_field->nr = equ->nw;
  mag_field->nz = equ->nh;
  mag_field->b0r0 = equ->rcenter * equ->bcenter;
  mag_field->r = malloc(mag_field->nr * sizeof(double));
  for (int i=0; i<mag_field->nr; i++)
  {
    mag_field->r[i] = equ->r[i];
  }
  mag_field->z = malloc(mag_field->nz * sizeof(double));
  for (int j=0; j<mag_field->nz; j++)
  {
    mag_field->z[j] = equ->z[j];
  }

  mag_field->Brz = allocate_3d_array(mag_field->nr, mag_field->nz, 2);
  double ***grad_psi_tmp = allocate_3d_array(mag_field->nr, mag_field->nz, 2);

  //calculate dpsi/dr and dpsi/dz at each point.
  selected_fun(mag_field->nr, mag_field->r, mag_field->nz, mag_field->z, equ->psi, grad_psi_tmp);
  
  for (int i=0; i<mag_field->nr; i++)
  {
    for (int j=0; j<mag_field->nz; j++)
    {
      if (mag_field->r[i] < MIN_R) //avoid 0 for division
      {
        mag_field->Brz[i][j][0] = -grad_psi_tmp[i][j][1]/( MIN_R);
        mag_field->Brz[i][j][1] =  grad_psi_tmp[i][j][0]/( MIN_R);
      }
      else
      {
        mag_field->Brz[i][j][0] = -grad_psi_tmp[i][j][1]/( mag_field->r[i]);
        mag_field->Brz[i][j][1] =  grad_psi_tmp[i][j][0]/( mag_field->r[i]);
      }
      
    }
  }
  //allocate and calculate dBrzdx, dBrzdx, d2Brzdxdy at each point by 2nd central difference method
  //ToDO in the future can be updated to high order accuracy.

  mag_field->dBrzdx = allocate_3d_array(mag_field->nr, mag_field->nz, 2);
  mag_field->dBrzdy = allocate_3d_array(mag_field->nr, mag_field->nz, 2);
  mag_field->d2Brzdxdy = allocate_3d_array(mag_field->nr, mag_field->nz, 2);
 
  int nr = mag_field->nr;
  int nz = mag_field->nz;
  for(int i=0;i<nr;i++)
    {
      for(int j=0;j<nz;j++)
      {
        int iyp=min(nz-1,j+1);
        int iym=max(0,j-1);
        int ixp=min(nr-1,i+1);
        int ixm=max(0,i-1);

        for (int k=0; k<2; k++)
        {
        mag_field->dBrzdx[i][j][k] = (mag_field->Brz[ixp][j][k] - mag_field->Brz[ixm][j][k])/(mag_field->r[ixp]-mag_field->r[ixm]);
        mag_field->dBrzdy[i][j][k] = (mag_field->Brz[i][iyp][k] - mag_field->Brz[i][iym][k])/(mag_field->z[iyp]-mag_field->z[iym]);
        mag_field->d2Brzdxdy [i][j][k] = (mag_field->Brz[ixp][iyp][k] - mag_field->Brz[ixm][iyp][k]
                                          -mag_field->Brz[ixp][iym][k] + mag_field->Brz[ixm][iym][k])
                                        /((mag_field->r[ixp]-mag_field->r[ixm])*(mag_field->z[iyp]-mag_field->z[iym]));
        }
      }
    }
    free_3d_array(grad_psi_tmp);
  }

void get_bt_torsys(MagFieldTorSys *mag_field, const double r0, double *bt)
{
  if (r0 < MIN_R)
  {
    *bt = mag_field->b0r0/MIN_R;
  }
  else
  {
    *bt = mag_field->b0r0/r0;
  }
}

void get_brz_torsys(MagFieldTorSys *mag_field, const double r, const double z, const char *method, double *br, double *bz)
{
  intpl_2d_fun selected_fun = get_intpl_2d_method(method);
  if (selected_fun == NULL)
  {
    printf("ERROE: Method: '%s' not found!\n", method);
    exit(1);
  }
  if (strcmp("bilnear", method)==0)
  {
    selected_fun(r, z, mag_field->nr, mag_field->r, mag_field->nz, mag_field->z, 
                mag_field->Brz, br, bz, NULL, NULL, NULL);
  }
  else{
    printf("Currently, do not support other methods!\n");
  }
  //need to do
  //else selected_fun(r, z, mag_field->nr, mag_field->r, mag_field->nz, mag_field->z, 
  //             mag_field->Brz, br, bz, mag_field);
} 


void write_mag_field_torsys(MagFieldTorSys *mag_field)
{
    write_brz_torsys(mag_field);
    write_bphi_torsys(mag_field);
    printf("Finish write magnetic field.\n");
    return;
}

void write_brz_torsys(MagFieldTorSys *mag_field)
{
  const char* filename;
  filename = "Brz_in_R";
  FILE* file1 = fopen(filename, "w");
  for (int i=0; i<mag_field->nr;i++)
  {
    for (int j = 0; j<mag_field->nz;j++)
    {
      fprintf(file1, "%lf  %lf  %lf\n", \
      mag_field->r[i], mag_field->z[j], mag_field->Brz[i][j][0]);

    }
  }
  fclose(file1);
  printf("write the Brz field in R direction in %s\n", filename);
//********************************************************************
  filename = "Brz_in_Z";
  FILE* file2 = fopen(filename, "w");
  for (int i=0; i<mag_field->nr;i++)
  {
    for (int j = 0; j<mag_field->nz;j++)
    {
      fprintf(file2, "%lf  %lf  %lf\n", \
      mag_field->r[i], mag_field->z[j], mag_field->Brz[i][j][1]);

    }
  }
  fclose(file2);
  printf("write the Brz field in Z direction in %s\n", filename);
//********************************************************************
  filename = "Bpol";
  FILE* file3 = fopen(filename, "w");
  double tmp;
  for (int i=0; i<mag_field->nr;i++)
  {
    for (int j = 0; j<mag_field->nz;j++)
    {
      
      tmp = mag_field->Brz[i][j][0] * mag_field->Brz[i][j][0] 
           +mag_field->Brz[i][j][1] * mag_field->Brz[i][j][1];
      tmp = sqrt(tmp);
      fprintf(file3, "%lf  %lf  %lf\n", \
      mag_field->r[i], mag_field->z[j], tmp);
    }
  }
  fclose(file3);
  printf("write the value of Bpol in %s\n", filename);
}

void write_bphi_torsys(MagFieldTorSys *mag_field)
{
  const char* filename;
  filename = "Btor";
  FILE* file = fopen(filename, "w");
  double tmp;
  for (int i=0; i<mag_field->nr;i++)
  {
    for (int j = 0; j<mag_field->nz;j++)
    {
      get_bt_torsys(mag_field,mag_field->r[i], &tmp);
      fprintf(file, "%lf  %lf  %lf\n", \
      mag_field->r[i], mag_field->z[j], tmp);
    }
  }
  fclose(file);
  printf("write the value of Btor in %s\n", filename);
}

