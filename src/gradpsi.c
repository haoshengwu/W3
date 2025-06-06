#include "gradpsi.h"
#include "mathbase.h"

GradPsiStr* init_grad_psi(void)
{
    GradPsiStr* gradpsi=malloc(sizeof(GradPsiStr));
    gradpsi->nr=-1;
    gradpsi->nz=-1;
    gradpsi->r=NULL;
    gradpsi->z=NULL;
    gradpsi->psi=NULL;
    gradpsi->gradpsi_r=NULL;
    gradpsi->gradpsi_z=NULL;
    gradpsi->dgradpsi_r_dr=NULL;
    gradpsi->dgradpsi_r_dz=NULL;
    gradpsi->dgradpsi_z_dr=NULL;
    gradpsi->dgradpsi_z_dz=NULL;
    gradpsi->d2gradpsi_r_drdz=NULL;
    gradpsi->d2gradpsi_z_drdz=NULL;
    return gradpsi;
}

void free_grad_psi(GradPsiStr *gradpsi)
{
    if(!gradpsi->r) free(gradpsi->r);
    if(!gradpsi->z) free(gradpsi->z);
    if(!gradpsi->psi) free_2d_array(gradpsi->psi);
    if(!gradpsi->gradpsi_r) free_2d_array(gradpsi->gradpsi_r);
    if(!gradpsi->gradpsi_z) free_2d_array(gradpsi->gradpsi_z);
    if(!gradpsi->dgradpsi_r_dr) free_2d_array(gradpsi->dgradpsi_r_dr);
    if(!gradpsi->dgradpsi_r_dz) free_2d_array(gradpsi->dgradpsi_r_dz);
    if(!gradpsi->dgradpsi_z_dr) free_2d_array(gradpsi->dgradpsi_z_dr);
    if(!gradpsi->dgradpsi_z_dz) free_2d_array(gradpsi->dgradpsi_z_dz);
    if(!gradpsi->d2gradpsi_r_drdz) free_2d_array(gradpsi->d2gradpsi_r_drdz);
    if(!gradpsi->d2gradpsi_z_drdz) free_2d_array(gradpsi->d2gradpsi_z_drdz);
}

void calc_grad_psi(Equilibrium *equ, GradPsiStr *gradpsi, diff2d_eval_fun func)
{
// func(nx, r[], ny, z[], f[nx][ny], dfdx[nx][ny], dfdy[nx][ny])
  int nr = equ->nw;
  int nz = equ->nh;
  gradpsi->nr = nr;
  gradpsi->nz = nz;
  gradpsi->psi=allocate_2d_array(nr,nz);
  gradpsi->r = malloc(gradpsi->nr*sizeof(double));
  gradpsi->z = malloc(gradpsi->nz*sizeof(double));

  for(int i=0; i<nr; i++)
  {
    gradpsi->r[i]=equ->r[i];
  }

  for(int j=0; j<nz; j++)
  {
    gradpsi->z[j]=equ->z[j];
  }

  for(int i=0; i<nr; i++)
  {
    for (int j=0; j<nz; j++)
    {
      gradpsi->psi[i][j]=equ->psi[i][j];
    }
  }

  gradpsi->gradpsi_r=allocate_2d_array(nr,nz);
  gradpsi->gradpsi_z=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_r_dr=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_r_dz=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_z_dr=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_z_dz=allocate_2d_array(nr,nz);
  gradpsi->d2gradpsi_r_drdz=allocate_2d_array(nr,nz);
  gradpsi->d2gradpsi_z_drdz=allocate_2d_array(nr,nz);
  func(nr, gradpsi->r, nz, gradpsi->z, gradpsi->psi, gradpsi->gradpsi_r, gradpsi->gradpsi_z);
  func(nr, gradpsi->r, nz, gradpsi->z,gradpsi->gradpsi_r, gradpsi->dgradpsi_r_dr, gradpsi->dgradpsi_r_dz);
  func(nr, gradpsi->r, nz, gradpsi->z,gradpsi->gradpsi_z, gradpsi->dgradpsi_z_dr, gradpsi->dgradpsi_z_dz);

  for (int i = 0; i < nr; i++) 
  {
    for (int j = 0; j < nz; j++) 
    {
      int ixp = min(nr-1, i+1);
      int ixm = max(0,     i-1);
      int iyp = min(nz-1, j+1);
      int iym = max(0,     j-1);

      double dr = gradpsi->r[ixp] - gradpsi->r[ixm];
      double dz = gradpsi->z[iyp] - gradpsi->z[iym];
      double scale = 1.0 / (dr * dz);

      gradpsi->d2gradpsi_r_drdz[i][j] = (
          gradpsi->gradpsi_r[ixp][iyp] - gradpsi->gradpsi_r[ixm][iyp]
        - gradpsi->gradpsi_r[ixp][iym] + gradpsi->gradpsi_r[ixm][iym]
      ) * scale;

      gradpsi->d2gradpsi_z_drdz[i][j] = (
        gradpsi->gradpsi_z[ixp][iyp] - gradpsi->gradpsi_z[ixm][iyp]
        - gradpsi->gradpsi_z[ixp][iym] + gradpsi->gradpsi_z[ixm][iym]
      ) * scale;
    }
  }
}

