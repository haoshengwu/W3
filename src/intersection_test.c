#include "intersection_test.h"
#include "curve.h"
#include "config.h"


PlateCells* generate_platecells(ThreeDimGrid* grid, int idx_zone)
{
  if(!grid)
  {
    fprintf(stderr, "Error: Null inputs for generate_platecells.\n");
    exit(EXIT_FAILURE);
  }
  PlateCells* platecells = malloc(sizeof(PlateCells));
  
  if (!platecells) {
    fprintf(stderr, "Error: failed to allocate memory for platecells.\n");
    exit(EXIT_FAILURE);
  }

  int np=grid->npol;
  int nr=grid->nrad;
  int nt=grid->ntor;

  int npcell=np-1;
  int nrcell=nr-1;
  int ntcell=nt-1;

  if(np<2||nr<2||nt<2)
  {
    fprintf(stderr, "Error: Invaild size of 3D grid in generate_platecells.\n");
    exit(EXIT_FAILURE);
  }

  platecells->npcell=npcell;
  platecells->nrcell=nrcell;
  platecells->ntcell=ntcell;

  platecells->idx_zone=idx_zone;

  platecells->is_plates=allocate_3d_array_int(nt-1,np-1,nr-1);

  #ifdef DEBUG
    for(int it=0;it<ntcell;it++)
    {
      for(int ip=0;ip<npcell;ip++)
      {
        for(int ir=0;ir<nrcell;ir++)
        {
          if(platecells->is_plates[it][ip][ir]!=0)
          {
            fprintf(stderr, "Error: Initialization value should be ZERO.\n");
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  #endif
  return platecells;
}

void free_platecells(PlateCells* platecells)
{
  if(platecells==NULL)
  {
    return;
  }
  free_3d_array_int(platecells->is_plates);
  platecells->is_plates=NULL;
  free(platecells);
  #ifdef DEBUG
    printf("Successfully free platecells.\n");
  #endif
}

void set_is_plate_cell(PlateCells* platecells, int ipcell, int ircell, int itcell)
{
  if (!platecells) {
    fprintf(stderr, "Error: Null platecells pointer set_is_plate_cell.\n");
    exit(EXIT_FAILURE);
  }
  if(ipcell<0 || ipcell>platecells->npcell-1)
  {
    fprintf(stderr, "Error: poloidal cell index is out of range in set_is_plate_cell.\n");
    exit(EXIT_FAILURE);
  }
  
  if(ircell<0 || ircell>platecells->nrcell-1)
  {
    fprintf(stderr, "Error: radial cell index is out of range in set_is_plate_cell.\n");
    exit(EXIT_FAILURE);
  }
  
  if(itcell<0 || itcell>platecells->ntcell-1)
  {
    fprintf(stderr, "Error: toroidal cell index is out of range in set_is_plate_cell.\n");
    exit(EXIT_FAILURE);
  }
  platecells->is_plates[itcell][ipcell][ircell]=1;
}

int get_plate_cell(PlateCells* platecells, int ipcell, int ircell, int itcell)
{
  if (!platecells) {
    fprintf(stderr, "Error: Null platecells pointer in get_plate_cell.\n");
    exit(EXIT_FAILURE);
  }
  if(ipcell<0 || ipcell>platecells->npcell-1)
  {
    fprintf(stderr, "Error: poloidal cell index is out of range in get_plate_cell.\n");
    exit(EXIT_FAILURE);
  }
  
  if(ircell<0 || ircell>platecells->nrcell-1)
  {
    fprintf(stderr, "Error: radial cell index is out of range in get_plate_cell.\n");
    exit(EXIT_FAILURE);
  }
  
  if(itcell<0 || itcell>platecells->ntcell-1)
  {
    fprintf(stderr, "Error: toroidal cell index is out of range in get_plate_cell.\n");
    exit(EXIT_FAILURE);
  }

  return(platecells->is_plates[itcell][ipcell][ircell]);
}


void update_platecells_radneu(PlateCells* platecells, int n_rad_neu)
{
  if(!platecells)
  {
    fprintf(stderr, "Error: Null inputs for update_platecells_radneu.\n");
    exit(EXIT_FAILURE);
  }

  if(n_rad_neu<0 || n_rad_neu>platecells->nrcell)
  {
    fprintf(stderr, "Error: Invalid radial neutral region index for update_platecells_radneu.\n");
    exit(EXIT_FAILURE);
  }

  int npcell=platecells->npcell;
  int nrcell=platecells->nrcell;
  int ntcell=platecells->ntcell;

  int idx_rad_cell=nrcell-n_rad_neu;

  for(int it=0;it<ntcell;it++)
  {
    for(int ip=0;ip<npcell;ip++)
    {
      for(int ir=idx_rad_cell; ir<nrcell; ir++)
      {
        set_is_plate_cell(platecells,ip,ir,it);
      }
    }
  }
  #ifdef DEBUG
  printf("Finish update the cells within/behind the neutral domain.\n");
  #endif
}

void update_platecells_targets(PlateCells* platecells, ThreeDimGrid* grid)
{
  //check input
  if(!platecells || !grid)
  {
    fprintf(stderr, "Error: Null inputs for update_platecells_targets.\n");
    exit(EXIT_FAILURE);
  }

  //check size
  int np=grid->npol;
  int nr=grid->nrad;
  int nt=grid->ntor;

  int npcell=platecells->npcell;
  int nrcell=platecells->nrcell;
  int ntcell=platecells->ntcell;

  if(npcell!=(np-1) || nrcell!=(nr-1) || ntcell!=(nt-1))
  {
    fprintf(stderr, "Error: The sizes of platecells and 3D grid are not match in update_platecells_targets.\n");
    exit(EXIT_FAILURE);
  }

  //The first and last cells are definitely is_place_cell
  for(int itcell=0;itcell<ntcell;itcell++)
  {
    for(int ircell=0;ircell<nrcell;ircell++)
    {
      set_is_plate_cell(platecells,0,ircell,itcell);
      set_is_plate_cell(platecells,npcell-1,ircell,itcell);
    }
  }
  //Filling the cells within/behind outer target plate
  for(int ipcell=1;ipcell<=ntcell;ipcell++)
  {
    for(int itcell=0;itcell<ntcell-ipcell+1;itcell++)
    {
      for(int ircell=0;ircell<nrcell;ircell++)
      {
        set_is_plate_cell(platecells,ipcell,ircell,itcell);
      }
    }
  }
  //Filling the cells within/behind inner target plate
  //i is the i-th(0-based) element from the end
  for(int i=1;i<=ntcell;i++)
  {
    for(int itcell=i-1;itcell<ntcell;itcell++)
    {
      for(int ircell=0;ircell<nrcell;ircell++)
      {
        set_is_plate_cell(platecells,npcell-1-i,ircell,itcell);
      }
    }
  }
  #ifdef DEBUG
  printf("Finish update the cells within/behind the inner&outer targets.\n");
  #endif
}

static void centroid4(double x1, double y1,
                      double x2, double y2,
                      double x3, double y3,
                      double x4, double y4,
                      double *cx, double *cy) 
{
    double A = 0.0, Cx = 0.0, Cy = 0.0;
    double cross;

    // Edge (x1,y1) → (x2,y2)
    cross = x1 * y2 - x2 * y1;
    A  += cross;
    Cx += (x1 + x2) * cross;
    Cy += (y1 + y2) * cross;

    // Edge (x2,y2) → (x3,y3)
    cross = x2 * y3 - x3 * y2;
    A  += cross;
    Cx += (x2 + x3) * cross;
    Cy += (y2 + y3) * cross;

    // Edge (x3,y3) → (x4,y4)
    cross = x3 * y4 - x4 * y3;
    A  += cross;
    Cx += (x3 + x4) * cross;
    Cy += (y3 + y4) * cross;

    // Edge (x4,y4) → (x1,y1)
    cross = x4 * y1 - x1 * y4;
    A  += cross;
    Cx += (x4 + x1) * cross;
    Cy += (y4 + y1) * cross;

    // Polygon area
    A *= 0.5;

    if (fabs(A) < EPSILON_12) 
    {
      *cx = 0.25*(x1+x2+x3+x4);
      *cy = 0.25*(y1+y2+y3+y4);
      return;
    }
    // Centroid coordinates
    *cx = Cx / (6.0 * A);
    *cy = Cy / (6.0 * A);
}

//For a 3D grid, at the it-th(0-base) phi-slice, 
// ipcell and ircell are the indexs of a cell.
// This function calculate the centroid4 of the cell
static void centroid4_of_cell(ThreeDimGrid* grid,
                              int ipcell, int ircell,
                              int it,
                              double *cx, double *cy)
{
  double x1=get_r_3Dgrid(grid, ipcell, ircell, it);
  double y1=get_z_3Dgrid(grid, ipcell, ircell, it);

  double x2=get_r_3Dgrid(grid, ipcell, ircell+1, it);
  double y2=get_z_3Dgrid(grid, ipcell, ircell+1, it);

  double x3=get_r_3Dgrid(grid, ipcell+1, ircell+1, it);
  double y3=get_z_3Dgrid(grid, ipcell+1, ircell+1, it);

  double x4=get_r_3Dgrid(grid, ipcell+1, ircell, it);
  double y4=get_z_3Dgrid(grid, ipcell+1, ircell, it);

  centroid4(x1,y1,x2,y2,x3,y3,x4,y4,cx,cy);

}
void update_platecells_limiter(PlateCells* platecells, ThreeDimGrid* grid, DLListNode* limiter_head, 
                               const double start_phi, const double end_phi)
{

  //check input
  if(!platecells || !grid)
  {
    fprintf(stderr, "Error: Null inputs for update_platecells_limiter.\n");
    exit(EXIT_FAILURE);
  }

  //check size
  int np=grid->npol;
  int nr=grid->nrad;
  int nt=grid->ntor;

  int npcell=platecells->npcell;
  int nrcell=platecells->nrcell;
  int ntcell=platecells->ntcell;

  if(npcell!=(np-1) || nrcell!=(nr-1) || ntcell!=(nt-1))
  {
    fprintf(stderr, "Error: The sizes of platecells and 3D grid are not match in update_platecells_limiter.\n");
    exit(EXIT_FAILURE);
  }

  //check limiter phi
  int idx_start_phi=-1;
  int idx_end_phi=-1;
  for(int it=0;it<nt;it++)
  {
    double phi=get_phi_3Dgrid(grid,0,0,it);
    if(fabs(start_phi-phi)<EPSILON_8)
    {
      idx_start_phi=it;
    }
    if(fabs(end_phi-phi)<EPSILON_8)
    {
      idx_end_phi=it;
    }
  }

  if(idx_start_phi==-1 || idx_end_phi==-1)
  {
    fprintf(stderr, "Error: The start/end phi of limiter is not match any phi in 3D grid in update_platecells_limiter.\n");
    exit(EXIT_FAILURE);
  }

  if(idx_start_phi>=idx_end_phi)
  {
    fprintf(stderr, "The phi range of limiter is wrong in update_platecells_limiter.\n");
    exit(EXIT_FAILURE);
  }

// The basic idea for filling plate cells is to check the centroid coordinates of each cell in phi-slice along the radial direction.
// If the centroid lies within the closed limiter structure, then all remaining radial cells are also considered inside the plates.
// In the toroidal direction, the logic is simpler: if the poloidal section of a cell in the phi-slice lies within the structure, 
// we assume the entire cell is contained as well.
  Curve* boundary[2];
  boundary[0]=convert_ddl_to_curve(limiter_head);
  boundary[1]=create_curve(0);
  add_last_point_curve(boundary[1],
                       get_curve_x(boundary[0],boundary[0]->n_point-1),
                       get_curve_y(boundary[0],boundary[0]->n_point-1));
  add_last_point_curve(boundary[1],
                       get_curve_x(boundary[0],0),
                       get_curve_y(boundary[0],0));
  #ifdef DEBUG
  write_curve("Limiter_bnd0",boundary[0]);
  write_curve("Limiter_bnd1",boundary[1]);
  #endif
  //for toroidal direction, currently it is simple
  //We only check the [idx_start_phi:idx_end_phi-1] phi-slice.
  for(int it=idx_start_phi;it<idx_end_phi-1;it++)
  {
    for(int ipcell=0;ipcell<npcell;ipcell++)
    {
      for(int ircell=0;ircell<nrcell;ircell++)
      {
        double cx;
        double cy;
        centroid4_of_cell(grid, ipcell, ircell, it, &cx, &cy);
        if(is_point_inside(cx,cy, boundary, 2))
        {
          set_is_plate_cell(platecells, ipcell, ircell, it);
          // #ifdef DEBUG
          // printf("The platecells ipcell ircell itcell is filled: %d %d %d.\n",ipcell, ircell, it);
          // #endif
        }
      }
    }
  }
  //for idx_end_phi phi-slice, corresponding to idx_end_phi-1 cell
  for(int ipcell=0;ipcell<npcell;ipcell++)
  {
    for(int ircell=0;ircell<nrcell;ircell++)
    {
      double cx;
      double cy;
      centroid4_of_cell(grid, ipcell, ircell, idx_end_phi, &cx, &cy);
      if(is_point_inside(cx,cy, boundary, 2))
      {
        set_is_plate_cell(platecells, ipcell, ircell, idx_end_phi-1);
        // #ifdef DEBUG
        // printf("The platecells ipcell ircell itcell is filled: %d %d %d.\n",ipcell, ircell, idx_end_phi-1);
        // #endif
      }
    }
  }
  free_curve(boundary[0]);
  free_curve(boundary[1]);
  #ifdef DEBUG
  printf("Finish update the cells within/behind the limiter plate.\n");
  #endif
}


//the start and end indexs along toroidal direction
//in this range, the plates are within/behind plates
//ONLY USED in this file (intersection_test.c)
typedef struct {
    int start;  // start index of a run of 1s (inclusive)
    int end;    // end index of a run of 1s (inclusive)
} CellRange;

static CellRange* find_cell_ranges(const int *arr, const int n, int *out_count) 
{
  if (!arr || n <= 0)
  {
    printf("WARNING: Invalid inputs for find_cell_ranges!\n");
    if (out_count) *out_count = 0;
    return NULL;
  }
  
  // First pass: count how many ranges
  int count = 0;
  for (int i = 0; i < n; i++)
  {
    if (arr[i] == 1 && (i == 0 || arr[i - 1] == 0)) count++;
  }
  if (out_count) *out_count = count;  
  if (count == 0) return NULL;

  // Allocate array to store all ranges
  CellRange *ranges = (CellRange*)malloc(sizeof(CellRange) * count);
  if (!ranges) 
  {
    fprintf(stderr, "Failed to allocated memory for the array of CellRange!\n");
    exit(EXIT_FAILURE);
  }

  // Second pass: fill each [start, end] range
  int idx = 0;
  for (int i = 0; i < n; i++) 
  {
    if (arr[i] == 1)
    {
      int start = i;
      while (i + 1 < n && arr[i + 1] == 1) i++;
      ranges[idx].start = start;
      ranges[idx].end   = i;
      idx++;
    }
  }
  return ranges;
}

void write_platecells(PlateCells* platecells, char* filename)
{
  if(!platecells || !filename)
  {
    fprintf(stderr, "Error: Null inputs for write_platecells.\n");
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

  int npcell=platecells->npcell;
  int nrcell=platecells->nrcell;
  int ntcell=platecells->ntcell;
  int idx_zone=platecells->idx_zone;

  if(npcell<2||nrcell<2||ntcell<2)
  {
    fprintf(stderr, "Error: Invalid size of platecells.\n");
    exit(EXIT_FAILURE);
  }

  const char *fmt="%-3d%6d%6d%6d";
  
  FILE *fp = fopen(filename, "w");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }

  //temporary array used to store the plate cells along toroidal direction
  int* array_tmp=malloc(ntcell*sizeof(int));
  if(!array_tmp)
  {
    fprintf(stderr, "Error: Failed to allocate memmory for array_tmp.\n");
    exit(EXIT_FAILURE);
  }

  for(int ircell=0;ircell<nrcell;ircell++)
  {
    for(int ipcell=0;ipcell<npcell;ipcell++)
    {
      for(int itcell=0;itcell<ntcell;itcell++)
      {
        array_tmp[itcell]=get_plate_cell(platecells,ipcell,ircell,itcell);
      }
      int nrange;
      CellRange* ranges=find_cell_ranges(array_tmp, ntcell, &nrange);
      if(nrange==0 && ranges==NULL)
      {
        continue;
      }
      //write one
      else if(nrange!=0 && ranges!=NULL)
      {
        fprintf(fp, fmt, idx_zone, ircell, ipcell,2*nrange);
        for(int i=0;i<nrange;i++)
        {
          fprintf(fp,"%6d%6d",ranges[i].start, ranges[i].end);
        }
        fprintf(fp,"\n");
        free(ranges);
      }
      else
      {
        fprintf(stderr, "UNEXPECTED ERROR: ranges is not consistent with nrange.\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  fclose(fp);
  free(array_tmp);
  #ifdef DEBUG
  printf("Successfully write platecells to %s.\n",filename);
  #endif
}