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

double get_r_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
    return g->points[idx_3Dgrid(g, ip, ir, it)].r;
}

double get_z_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
    return g->points[idx_3Dgrid(g, ip, ir, it)].z;
}

double get_phi_3Dgrid(const ThreeDimGrid* g, int ip, int ir, int it)
{
  return g->points[idx_3Dgrid(g, ip, ir, it)].phi;
}

void set_point_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double r, double z, double phi)
{
    GridPoint3D* p = &g->points[idx_3Dgrid(g, ip, ir, it)];
    p->r = r;
    p->z = z;
    p->phi = phi;
}

void set_r_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double r)
{
  GridPoint3D* p = &g->points[idx_3Dgrid(g, ip, ir, it)];
  p->r = r;
}

void set_z_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double z)
{
  GridPoint3D* p = &g->points[idx_3Dgrid(g, ip, ir, it)];
  p->z = z;
}

void set_phi_3Dgrid(ThreeDimGrid* g, int ip, int ir, int it, double phi)
{
  GridPoint3D* p = &g->points[idx_3Dgrid(g, ip, ir, it)];
  p->phi = phi;
}

void assign_2D_to_3D_tor_slice(const TwoDimGrid* grid2d, ThreeDimGrid* grid3d, int it, double phim)
{
  //Check 
  if(!grid2d || !grid3d || !grid2d->points || !grid3d->points)
  {
    fprintf(stderr, "Error: Invalid inputs for assign_2D_to_3d_tor_slice.\n");
    exit(EXIT_FAILURE);
  }
  if(it<0 || it>grid3d->ntor-1)
  {
    fprintf(stderr, "Error: Toroidal index it is out of range in assign_2D_to_3D_tor_slice.\n");
    exit(EXIT_FAILURE);
  }
  if(grid2d->npol!=grid3d->npol||grid2d->nrad!=grid3d->nrad)
  {
    fprintf(stderr, "Error: The 2D grid size is not consistent with 3D grid slice.\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0;i<grid3d->npol;i++)
  {
    for(int j=0;j<grid3d->nrad;j++)
    {
      double r = get_x_2Dgrid(grid2d,i,j);
      double z = get_y_2Dgrid(grid2d,i,j);
      if(fabs(r)<EPSILON_12 && fabs(z)<EPSILON_12)
      {
        printf("WARNING: The point is %.12f %.12f\n", r, z);
      }
      set_point_3Dgrid(grid3d,i,j,it,r,z,phim);
    }
  }
}

void generate_EMC3_3Dgrid_from_2Dgrid_tracing(const TwoDimGrid* grid2d, ThreeDimGrid* grid3d, 
                                              double phim, int nphi, double* phi,
                                              ode_function* func,ode_solver* solver)
{
  /*****************
  *  Check inputs  *
  ******************/
  if(!grid2d||!grid3d)
  {
    fprintf(stderr, "Error: Invalid inputs generate_EMC3_3Dgrid_from_2Dgrid_tracing.\n");
    exit(EXIT_FAILURE);
  }

  if(grid2d->npol!=grid3d->npol||grid2d->nrad!=grid3d->nrad)
  {
    fprintf(stderr, "Error: The poloidal and radial sizes of grid2D and grid3D are not identical.\n");
    exit(EXIT_FAILURE);
  }
  
  int np=grid2d->npol;
  int nr=grid2d->nrad;
  int nt=grid3d->ntor;

  if(nphi!=nt)
  {
    fprintf(stderr, "Error: The toroidal size of grid3D is not identical with nphi.\n");
    exit(EXIT_FAILURE);
  }

  if(func->ndim!=3)
  {
    fprintf(stderr, "Error: The 3D magnetic field is needed.\n");
    exit(EXIT_FAILURE);
  }
  
  if(grid3d->opt_direction!=GRID_3D_OPTIMIZE_RPT)
  {
    fprintf(stderr, "Error: EMC3 3D GRID MUST BE RAD-POL-TOR ORDER。\n");
    exit(EXIT_FAILURE);
  }
  int nfirst, nlast, idx_phim;
  calc_nfirst_nlast(phim, nphi, phi, &nfirst, &nlast);
  idx_phim=nlast;

  /***********************************
  *  Check magnetic field direction  *
  ***********************************/
  restore_3D_mag_direction(func);
  double pt[3], pt_next[3];
  //using the point ip=1&ir=1, because of ip=0 ir=0 can be X-poiint
  pt[0]=get_x_2Dgrid(grid2d,1,1);
  pt[1]=get_y_2Dgrid(grid2d,1,1);
  pt[2]=phim;
  double t_tmp=0.0;
  solver->next_step(solver->step_size, &t_tmp, pt, pt_next, solver->solver_data, func);
  if ((phi[0]-phim) * (pt_next[2] - pt[2]) < 0)
  {
    fprintf(stderr, "Error: direction phi1->phi2 is not same with magnetic field direction.\n");
    exit(EXIT_FAILURE);
  }
  #ifdef DEBUG
    printf("The direction of magnetic field is ok.\n");
  #endif
  /*****************************************
  *  assign grid2d to phi[idx_phim] plane  *
  ******************************************/
  assign_2D_to_3D_tor_slice(grid2d, grid3d, idx_phim, phim);
  
  #ifdef DEBUG
    printf("Assign grid2d to Toroidal %d (0-base) slice to grid3d\n",idx_phim);
  #endif

  /***********************************************************
  *  Generate the grid from phi[idx_phim-1] plane to phi[0]  *
  ************************************************************/
  TwoDimGrid* grid2d_tmp=create_2Dgrid_radial_major(np,nr);
  restore_3D_mag_direction(func);
  for(int i=idx_phim-1;i>-1;i--) //i is index!
  {
    generate_2Dgrid_tracing(grid2d, phim, grid2d_tmp, phi[i], func, solver);
    assign_2D_to_3D_tor_slice(grid2d_tmp, grid3d, i, phi[i]);
    #ifdef DEBUG
      printf("Assign grid2d to Toroidal %d (0-base) slice to grid3d\n",i);
    #endif
  }

  #ifdef DEBUG
    printf("Finish the grid3d from phi[%d] plane to phi[%d]\n",idx_phim-1,0);
  #endif
  /*************************************************************
  *  Generate the grid from phi[idx_phim+1] plane to phi[nt-1] *
  *************************************************************/
  reverse_3D_mag_direction(func);
  for(int i=idx_phim+1;i<nt;i++) //i is index!
  {
    generate_2Dgrid_tracing(grid2d, phim, grid2d_tmp, phi[i], func, solver);
    assign_2D_to_3D_tor_slice(grid2d_tmp, grid3d, i, phi[i]);
    #ifdef DEBUG
      printf("Assign grid2d to Toroidal %d (0-base) slice to grid3d\n",i);
    #endif
  }

  #ifdef DEBUG
    printf("Finish the grid3d from phi[%d] plane to phi[%d]\n",idx_phim+1,nt-1);
  #endif

  free_2Dgrid(grid2d_tmp);
}

void write_EMC3_3Dgrid_to_XYZ_CSYS(ThreeDimGrid* g, char* filename)
{
  if (!filename || !g) 
  {
    fprintf(stderr, "Error: NULL input to write_EMC3_3Dgrid_to_XYZ_CSYS.\n");
    exit(EXIT_FAILURE);
  }

  FILE *fp = fopen(filename, "w");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }
  int np=g->npol;
  int nr=g->nrad;
  int nt=g->ntor;
  fprintf(fp, "# %d %d %d %s\n",np, nr, nt,
         (g->opt_direction == GRID_3D_OPTIMIZE_RPT) ? "GRID_3D_OPTIMIZE_RPT" : "GRID_3D_OPTIMIZE_PRT");

  for(int it=0; it<nt; it++)
  {
    double phi=get_phi_3Dgrid(g,0,0,it);
    fprintf(fp, "# PHI %d %.12f\n",it, phi);
    for(int ip=0; ip<np; ip++)
    {
      for(int ir=0; ir<nr; ir++)
      {
        double r=get_r_3Dgrid(g,ip,ir,it);
        //The Z from RZPHI_CSYS is same with XYZ_CSYS.
        double z=get_z_3Dgrid(g,ip,ir,it);
        double x=r * cos(deg2rad(phi));
        double y=r * sin(deg2rad(phi));
        fprintf(fp, "%.12f %.12f %.12f\n", x, y, z);
      }
    }
  }
  fclose(fp);
}

void write_EMC3_3Dgrid_to_EMC3_format(ThreeDimGrid* g, char* filename)
{
  if (!filename || !g) 
  {
    fprintf(stderr, "Error: NULL input to write_EMC3_3Dgrid_to_EMC3_format.\n");
    exit(EXIT_FAILURE);
  }

  FILE *fp = fopen(filename, "w");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }
  int np=g->npol;
  int nr=g->nrad;
  int nt=g->ntor;
  fprintf(fp, "%10d%10d%10d\n",nr, np, nt);
  for(int it=0; it<nt; it++)
  {
    double phi=get_phi_3Dgrid(g,0,0,it);
    fprintf(fp, "%15.8f\n",phi);
    for(int ip=0; ip<np; ip++)
    {
      for(int ir=0; ir<nr; ir++)
      {
        double r=get_r_3Dgrid(g,ip,ir,it);
        fprintf(fp, "%.15e\n", r);
      }
    }
    for(int ip=0; ip<np; ip++)
    {
      for(int ir=0; ir<nr; ir++)
      {
        double z=get_z_3Dgrid(g,ip,ir,it);
        fprintf(fp, "%.15e\n", z);
      }
    }
  }
  fclose(fp);
  #ifdef DEBUG
  printf("Successfully write EMC3 format 3D GIRD: %s.\n",filename);
  #endif
}

ThreeDimGrid* load_EMC3_format_3Dgrid_from_file(char* filename)
{
  if (!filename) 
  {
    fprintf(stderr, "Error: NULL input to load_EMC3_format_3Dgrid_from_file.\n");
    exit(EXIT_FAILURE);
  }
 
  
  FILE *fp = fopen(filename, "r");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }
  

  int np;
  int nr;
  int nt;

  if(fscanf(fp, "%10d%10d%10d", &nr,&np,&nt)!=3)
  {
    fprintf(stderr, "Error: Cannot read the size of the 3D grid in %s.\n",filename);
    fclose(fp);
    exit(EXIT_FAILURE);
  }
  else
  {
    printf("Begin to read the 3D grid from %s.\n",filename);
    printf("The size of 3D grid is nr np nt %6d %6d %6d.\n", nr, np, nt);
  }

  ThreeDimGrid* grid3d=create_3Dgrid_radial_major(np,nr,nt);

  double phi;
  for(int it=0;it<nt;it++)
  {
    fscanf(fp, "%lf", &phi);
    #ifdef DEBUG
    printf("Read the toroidal index %6d at phi %12f\n",it,phi);
    #endif
    double r;
    double z;
    for(int ip=0; ip<np; ip++)
    {
      for(int ir=0; ir<nr; ir++)
      {
        fscanf(fp, "%lf", &r);
        set_r_3Dgrid(grid3d, ip,ir,it, r);
        set_phi_3Dgrid(grid3d, ip,ir,it, phi);
      }
    }
    for(int ip=0; ip<np; ip++)
    {
      for(int ir=0; ir<nr; ir++)
      {
        fscanf(fp, "%lf", &z);
        set_z_3Dgrid(grid3d, ip,ir,it, z);
      }
    }
  }
  if (fscanf(fp, " %lf", &phi) != EOF) 
  {
    fprintf(stderr, "Unexpected extra data in file %s.\n", filename);
    fclose(fp);
    free_3Dgrid(grid3d);
    exit(EXIT_FAILURE);
  }
  fclose(fp);
  printf("Successfully read 3d grid from %s.\n", filename);
  return grid3d;
}