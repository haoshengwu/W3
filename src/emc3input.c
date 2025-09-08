#include "emc3input.h"

#include <math.h>
#include "config.h"
void write_BFIELD_file_default(ThreeDimGrid* grid, MagFieldTorSys* mag_field, char* filename)
{
  if(!grid || !mag_field || !filename)
  {
    fprintf(stderr, "Error: Null inputs for write_BFIELD_file_default.\n");
    exit(EXIT_FAILURE);
  }

  int np=grid->npol;
  int nr=grid->nrad;
  int nt=grid->ntor;

  double r;
  double z;
  
  FILE *fp = fopen(filename, "w");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }

  if (np <= 0 || nr <= 0 || nt <= 0) 
  {
    fprintf(stderr, "Error: Non-positive grid sizes (np=%d, nr=%d, nt=%d).\n", np, nr, nt);
    exit(EXIT_FAILURE);
  }

  for(int it=0;it<nt;it++)
  {
    for(int ip=0;ip<np;ip++)
    {
      for(int ir=0;ir<nr;ir++)
      {
        r = get_r_3Dgrid(grid, ip, ir, it);
        z = get_z_3Dgrid(grid, ip, ir, it);
        double br;
        double bz;
        //Calculate br, bz
        // cubicherm2d2f(r, z, mag_field->nr, mag_field->r, mag_field->nz, mag_field->z, 
        //               mag_field->Brz, &br, &bz, mag_field->dBrzdx, mag_field->dBrzdy, mag_field->d2Brzdxdy);
        bilenar2d2f(r, z, mag_field->nr, mag_field->r, mag_field->nz, mag_field->z, 
                      mag_field->Brz, &br, &bz, mag_field->dBrzdx, mag_field->dBrzdy, mag_field->d2Brzdxdy);
        //Calculate bt
        double bt = (fabs(r) > EPSILON_12) ? (mag_field->b0r0 / r) : 0.0;
        double bfield=sqrt(br*br+bz*bz+bt*bt);
        if (fprintf(fp, "%.15e\n", bfield) < 0) 
        {
          fprintf(stderr, "Error: write failed for file \"%s\".\n", filename);
          fclose(fp);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  fclose(fp);
  #ifdef DEBUG
  printf("Successfully write BFIELD of a 3D grid in : %s.\n",filename);
  #endif
}

void write_PLATE_MAG_file_test(ThreeDimGrid* grid, int n_neu_rad_cells, int pol_dir, int idx_zone, char* filename)
{
  if(!grid || !filename)
  {
    fprintf(stderr, "Error: Null inputs for write_PLATE_MAG_file_test.\n");
    exit(EXIT_FAILURE);
  }

  if(idx_zone<1)
  {
    printf("WARNING: Zone index 0 means CORE.\n");
    printf("WARNING: ARE YOU SURE?\n");
  }

  if(n_neu_rad_cells<0)
  {
    fprintf(stderr, "Error: Neutral region in the radial direction shoud >=0.\n");
    exit(EXIT_FAILURE);
  }

  //The numbers of poloidal, radial and toroidal points
  int np = grid->npol;
  int nr = grid->nrad;
  int nt = grid->ntor;

  //The numbers of poloidal, radial and toroidal cells
  int npcell=np-1;
  int nrcell=nr-1;
  int ntcell=nt-1;
 

  // Define plasma and neutral radial ranges
  int plasma_rad_range[2];
  int neu_rad_range[2];

  if (n_neu_rad_cells == 0) {
    // Special case: no neutral cells
    // → Plasma region covers all cells [0, nrcell-1]
    // → Neutral region is empty (set neu[0] > neu[1] to indicate an empty interval)
    plasma_rad_range[0] = 0;
    plasma_rad_range[1] = nrcell - 1;
    neu_rad_range[0]    = 1;   // start > end → empty range
    neu_rad_range[1]    = 0;
  } else {
    // General case: split into plasma and neutral regions
    // Plasma region: [0, nrcell - n_neutral_cells - 1]
    // Neutral region: [plasma_end + 1, nrcell - 1]
    plasma_rad_range[0] = 0;
    plasma_rad_range[1] = nrcell - n_neu_rad_cells - 1;
    neu_rad_range[0]    = plasma_rad_range[1] + 1;
    neu_rad_range[1]    = nrcell - 1;
  }

  FILE *fp = fopen(filename, "w");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }

  //EMC3 PLATE_MAG format: 
  /*
    For each line in the file: 
    * Index of zone (integer), cell index in radial (integer) and poloidal direction (integer). 
    * Number of values to read (integer). 
    * For each two value to read: 
    * Start (integer) and stop cell index in toroidal direction (integer). 
  */
  const char *fmt="%-3d%6d%6d%6d%6d%6d\n";

  if(pol_dir==0) //from inner to outer
  {
    for(int ir=0; ir<nrcell; ir++)
    {
      for(int ip=0; ip<npcell; ip++)
      {
        if(ir>=neu_rad_range[0] && ir<=neu_rad_range[1])  //Neutral domain
        {
          fprintf(fp,fmt,idx_zone,ir,ip,2,0,ntcell-1);
        }
        else if (ir>=plasma_rad_range[0] && ir<=plasma_rad_range[1]) //plasma domain
        {
          if(ip==0 || ip==ntcell-1)
          {
            fprintf(fp,fmt,idx_zone,ir,ip,2,0,ntcell-1);
          }
          if(ip>=1&&ip<=ntcell-1+1)
          { 
            //POLOIDAL [0:ntcell-1] cells, 0-base idx
            fprintf(fp,fmt,idx_zone,ir,ip,2,ip-1,ntcell-1); 
          } 
          else if (ip>=npcell-ntcell-1&&ip<=npcell-1-1)
          {
            //POLOIDAL [npcell-ntcell:npcell-1] cells, 0-base idx
            fprintf(fp,fmt,idx_zone,ir,ip,2,0,ip-(npcell-1-ntcell)); 
          }
        }
        else
        {
          fclose(fp);
          fprintf(stderr, "Unexpected Error: radial cell index %d out of range %d %d.\n", 
                  ir, plasma_rad_range[0], neu_rad_range[1]);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  else if(pol_dir==1) //from outer to inner
  {
    for(int ir=0; ir<nrcell; ir++)
    {
      for(int ip=0; ip<npcell; ip++)
      {
        if(ir>=neu_rad_range[0] && ir<=neu_rad_range[1]) //Neutral domain
        {
          fprintf(fp,fmt,idx_zone,ir,ip,2,0,ntcell-1);
        }
        else if (ir>=plasma_rad_range[0] && ir<=plasma_rad_range[1]) //plasma domain
        {
          if(ip==0 || ip==ntcell-1)
          {
            fprintf(fp,fmt,idx_zone,ir,ip,2,0,ntcell-1);
          }

          if(ip>=1&&ip<=ntcell-1+1)
          { 
            //POLOIDAL [0:ntcell-1] cells, 0-base idx
            fprintf(fp,fmt,idx_zone,ir,ip,2,0,ntcell-1-ip-1);
          } 
          else if (ip>=npcell-ntcell-1 && ip<=npcell-1-2)
          {
            //POLOIDAL [npcell-ntcell:npcell-1] cells, 0-base idx
            fprintf(fp,fmt,idx_zone,ir,ip,2,npcell-1-ip-1,ntcell-1);
          }
        }
        else
        {
          fclose(fp);
          fprintf(stderr, "Unexpected Error: radial cell index %d out of range %d %d.\n", 
                  ir, plasma_rad_range[0], neu_rad_range[1]);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  else
  {
    fclose(fp);
    fprintf(stderr, "Error: UNKNOWN Poloidal Direction: %d\n",pol_dir);
    exit(EXIT_FAILURE);
  }


}


void write_axis_sys_surface_default(ThreeDimGrid* grid, DLListNode* surface, bool reverse, char* filename)
{
  if(!grid || !filename || !surface)
  {
    fprintf(stderr, "Error: Null inputs for write_axis_sys_surface_default.\n");
    exit(EXIT_FAILURE);
  }
  int nt=grid->ntor;

  FILE *fp = fopen(filename, "w");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }


  DLListNode* surf_dll=copy_DLList(surface);
  if(reverse)
  {
    reverse_DLList(&surf_dll);
  }

  Curve* surf_curve=convert_ddl_to_curve(surf_dll);

  //write the title
  fprintf(fp,"#%s %f %f\n",filename, get_phi_3Dgrid(grid,0,0,0),get_phi_3Dgrid(grid,0,0,nt-1));
  int n_point=surf_curve->n_point;
  fprintf(fp,"%6d%6d%6d%12.6f%12.6f\n",nt, n_point, 1, 0.0, 0.0);
  for(int it=0;it<nt;it++)
  {
    fprintf(fp,"%12.6f\n", get_phi_3Dgrid(grid,0,0,it));
    for(int i=0;i<n_point;i++)
    {
      double r=surf_curve->points[i].x * 100; //FROM M to CM
      double z=surf_curve->points[i].y * 100; //FROM M to CM
      fprintf(fp,"%18.10f %18.10f\n",r,z);
    }
  }

  fclose(fp);
  #ifdef DEBUG
  printf("Successfully write 360 degree Axisymmetric Surface in: %s.\n",filename);
  #endif



  free_DLList(surf_dll);
  free_curve(surf_curve);
}

