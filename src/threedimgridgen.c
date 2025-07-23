#include "config.h"

#include "threedimgridgen.h"



ThreeDimGrid* create_3Dgrid_poloidal_major(int npol, int nrad, int ntor) 
{
    return create_3Dgrid_optimized_for(npol, nrad, ntor, GRID_3D_OPTIMIZE_PRT);
}

ThreeDimGrid* create_3Dgrid_radial_major(int npol, int nrad, int ntor) 
{
    return create_3Dgrid_optimized_for(npol, nrad, ntor, GRID_3D_OPTIMIZE_RPT);
}

ThreeDimGrid* create_3Dgrid_optimized_for(int npol, int nrad, int ntor, Grid3DOptimization opt) 
{
    // Validate input
    if (npol <= 0 || nrad <= 0 || ntor <= 0) 
    {
        fprintf(stderr, "Error: Invalid grid dimensions (%d, %d, %d)\n", npol, nrad, ntor);
        return NULL;
    }
    
    // Allocate grid structure
    ThreeDimGrid* grid = (ThreeDimGrid*)calloc(1, sizeof(ThreeDimGrid));
    if (!grid) 
    {
        fprintf(stderr, "Error: Failed to allocate grid structure\n");
        return NULL;
    }
    
    // Initialize grid parameters
    grid->npol = npol;
    grid->nrad = nrad;
    grid->ntor = ntor;
    
    int margin_pol = DEFAULT_POL_MARGIN;
    int margin_rad = DEFAULT_RAD_MARGIN;
    int margin_tor = DEFAULT_TOR_MARGIN;

    // Initialize with no extra capacity and zero offsets
    grid->cap_npol = npol+2*margin_pol;
    grid->cap_nrad = nrad+2*margin_rad;
    grid->cap_ntor = ntor+2*margin_tor;
    
    grid->offset_pol = margin_pol;
    grid->offset_rad = margin_rad;
    grid->offset_tor = margin_tor;
    
    grid->opt_direction = opt;
    
    // Calculate total points with overflow check
    size_t total_points = (size_t)grid->cap_npol * grid->cap_nrad * grid->cap_ntor;
    // Correct overflow check
    if (total_points / grid->cap_npol / grid->cap_nrad != grid->cap_ntor) 
    {
        fprintf(stderr, "Error: Grid dimensions too large (overflow)\n");
        free(grid);
        exit(EXIT_FAILURE);
    }
    
    // Allocate points array with alignment
    size_t alloc_size = total_points * sizeof(GridPoint3D);
    
    // C11 aligned_alloc requires size to be multiple of alignment
    size_t aligned_size = ((alloc_size + GRID_ALIGNMENT - 1) / GRID_ALIGNMENT) * GRID_ALIGNMENT;
    
    grid->points = (GridPoint3D*)aligned_alloc(GRID_ALIGNMENT, aligned_size);
    if (!grid->points) 
    {
        fprintf(stderr, "Error: Failed to allocate grid points (%zu bytes)\n", alloc_size);
        free(grid);
        exit(EXIT_FAILURE);
    }
    
    // Initialize all points to zero
    memset(grid->points, 0, aligned_size);  // 使用 aligned_size 而不是 alloc_size
    
    return grid;
}

void free_3Dgrid(ThreeDimGrid* grid) 
{
    if (!grid) 
    {
        return;  // Nothing to free
    }

    // Free points array if it exists
    if (grid->points) 
    {
        // Clear sensitive data if needed (optional)
        free(grid->points);
        grid->points = NULL;
    }

    // Free the grid structure
    free(grid);
}

static inline size_t idx_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
    if (g->opt_direction == GRID_3D_OPTIMIZE_PRT) 
    {
      // pol-rad-tor order (ip fastest, then ir, then it)
      return (size_t)(g->offset_tor + it) * g->cap_npol * g->cap_nrad
           + (size_t)(g->offset_rad + ir) * g->cap_npol
           + (size_t)(g->offset_pol + ip);
    } 
    else if(g->opt_direction == GRID_3D_OPTIMIZE_RPT)
    {
        // rad-pol-tor order (ir fastest, then ip, then it)
      return (size_t)(g->offset_tor + it) * g->cap_npol * g->cap_nrad
           + (size_t)(g->offset_pol + ip) * g->cap_nrad
           + (size_t)(g->offset_rad + ir);
    }
    else
    {
      fprintf(stderr, "Error: opt_direction in idx_3Dgrid.\n");
      exit(EXIT_FAILURE);
    }
}

static inline bool is_valid_index(const ThreeDimGrid* g, int ip, int ir, int it) 
{
    return (ip >= 0 && ip < g->npol 
         && ir >= 0 && ir < g->nrad 
         && it >= 0 && it < g->ntor);
}

GridPoint3D* get_point_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
    return &g->points[idx_3Dgrid(g, ip, ir, it)];
}

double get_x_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
    return g->points[idx_3Dgrid(g, ip, ir, it)].x;
}

double get_y_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
    return g->points[idx_3Dgrid(g, ip, ir, it)].y;
}

double get_z_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
  return g->points[idx_3Dgrid(g, ip, ir, it)].z;
}

void set_point_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double x, double y, double z)
{
    GridPoint3D* p = &g->points[idx_3Dgrid(g, ip, ir, it)];
    p->x = x;
    p->y = y;
    p->z = z;
}

void assign_2D_to_3d_tor_slice(const TwoDimGrid* grid2d, ThreeDimGrid* grid3d, int it)
{
  //Check 
  if(!grid2d || !grid3d || !grid2d->points || !grid3d->points)
  {
    fprintf(stderr, "Error: Invalid inputs for assign_2D_to_3d_tor_slice.\n");
    exit(EXIT_FAILURE);
  }
  if(it<0 || it>=grid3d->ntor-1)
  {
    fprintf(stderr, "Error: Toroidal index it is out of range.\n");
    exit(EXIT_FAILURE);
  }
  if(grid2d->npol!=grid3d->cap_npol||grid2d->cap_nrad!=grid3d->cap_nrad)
  {
    fprintf(stderr, "Error: The 2D grid size is not consistent with 3D grid slice.\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0;i<grid3d->npol;i++)
  {
    for(int j=0;j<grid3d->nrad;j++)
    {
      double x = get_x_2Dgrid(grid2d,i,j);
      double y = get_y_2Dgrid(grid2d,i,j);
      if(fabs(x)<EPSILON && fabs(y)<EPSILON)
      {
        printf("WARNING: The point is %.12f %.12f\n", x, y);
      }
      double z = get_z_3Dgrid(grid3d, i, j, it);
      set_point_3Dgrid(grid3d,i,j,it,x,y,z);
    }
  }
}
