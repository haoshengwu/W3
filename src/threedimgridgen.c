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


void expand_3Dgrid(ThreeDimGrid* g,
                   int add_pol_head, int add_pol_tail,
                   int add_rad_head, int add_rad_tail,
                   int add_tor_head, int add_tor_tail)
{
  if (!g) return;

  // Expansion values must be non-negative
  if (add_pol_head < 0 || add_pol_tail < 0 ||
      add_rad_head < 0 || add_rad_tail < 0 ||
      add_tor_head < 0 || add_tor_tail < 0)
  {
    fprintf(stderr, "expanded number should >= 0 for expand_3Dgrid\n");
    exit(EXIT_FAILURE);
  }

  // Compute new logical sizes
  int new_npol = g->npol + add_pol_head + add_pol_tail;
  int new_nrad = g->nrad + add_rad_head + add_rad_tail;
  int new_ntor = g->ntor + add_tor_head + add_tor_tail;

  // Compute new offsets (expanding at the head reduces offset index)
  int new_offset_pol = g->offset_pol - add_pol_head;
  int new_offset_rad = g->offset_rad - add_rad_head;
  int new_offset_tor = g->offset_tor - add_tor_head;

  // Boundary check: expansion must not exceed allocated capacity
  if (new_offset_pol < 0 || new_offset_rad < 0 || new_offset_tor < 0 ||
      (new_offset_pol + new_npol) > g->cap_npol ||
      (new_offset_rad + new_nrad) > g->cap_nrad ||
      (new_offset_tor + new_ntor) > g->cap_ntor)
  {
    fprintf(stderr,
      "expand_3Dgrid: expansion exceeds allocated memory (capacity: "
      "pol=%d, rad=%d, tor=%d). "
      "Requested offset_pol=%d, npol=%d, "
      "offset_rad=%d, nrad=%d, "
      "offset_tor=%d, ntor=%d\n",
      g->cap_npol, g->cap_nrad, g->cap_ntor,
      new_offset_pol, new_npol,
      new_offset_rad, new_nrad,
      new_offset_tor, new_ntor);
    fprintf(stderr, "In the future, dynamic expansion will be supported.\n");
    exit(EXIT_FAILURE);
  }

  // Apply updates to grid metadata
  g->offset_pol = new_offset_pol;
  g->offset_rad = new_offset_rad;
  g->offset_tor = new_offset_tor;

  g->npol = new_npol;
  g->nrad = new_nrad;
  g->ntor = new_ntor;

  // Note: underlying memory layout (points[]) depends on g->opt_direction.
  // No data move is needed as long as expansion stays within allocated capacity.
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


ThreeDimGrid* merge_3Dgrids_radial(const ThreeDimGrid* grid1, const ThreeDimGrid* grid2)
{
  if(!grid1 || !grid2)
  {
    fprintf(stderr, "Error: Invalid inputs for merge_3Dgrids_radial.\n");
    exit(EXIT_FAILURE);
  }

  if(grid1->npol!=grid2->npol || grid1->ntor!=grid2->ntor)
  {
    fprintf(stderr, "Error: The Pol or Tor size are not consistent for merge_3Dgrids_radial.\n");
    exit(EXIT_FAILURE);
  }
  //check the toroidal phi
  int ntor=grid1->ntor;
  for(int i=0;i<ntor;i++)
  {
    if (fabs(get_phi_3Dgrid(grid1,0,0,i)-fabs(get_phi_3Dgrid(grid2,0,0,i))>EPSILON_10))
    {
      fprintf(stderr, "In the toroidal direction, the phi are not exactly match.\n");
      exit(EXIT_FAILURE);
    }
  }
  //Check whether the last radial slice of grid1 is same with first slice of grid2.
  int npol=grid1->npol;
  int nrad1=grid1->nrad;
  int nrad2=grid2->nrad;
  int overlap=0;
  //check the frist points
  if(fabs(get_r_3Dgrid(grid1,0,nrad1-1,0)-get_r_3Dgrid(grid2,0,0,0))<EPSILON_10 &&
     fabs(get_z_3Dgrid(grid1,0,nrad1-1,0)-get_z_3Dgrid(grid2,0,0,0))<EPSILON_10)
  {
    overlap=1;
  }
  //the whole radial surface should be consistent with the first point.
  for(int k=0;k<ntor;k++)
  {
    for(int i=0;i<npol;i++)
    {
      int tmp=0;
      if(fabs(get_r_3Dgrid(grid1,i,nrad1-1,k)-get_r_3Dgrid(grid2,i,0,k))<EPSILON_10 &&
         fabs(get_z_3Dgrid(grid1,i,nrad1-1,k)-get_z_3Dgrid(grid2,i,0,k))<EPSILON_10)
      {
        tmp=1;
      }
      if(tmp!=overlap)
      {
        #ifdef DEBUG
        fprintf(stderr, " Position pol: %d, tor: %d.\n",i,k);
        fprintf(stderr, "grid1 r: %.15f z: %.15f\n",get_r_3Dgrid(grid1,i,nrad1-1,k), get_z_3Dgrid(grid1,i,nrad1-1,k));
        fprintf(stderr, "grid2 r: %.15f z: %.15f\n",get_r_3Dgrid(grid2,i,0,k), get_z_3Dgrid(grid2,i,0,k));
        #endif
        fprintf(stderr, "UNEXPECTED ERROR, the last radial slice of grid1 intersects with first slice of grid2?\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  ThreeDimGrid* new_grid=NULL;
  int nrad=nrad1+nrad2-overlap;
  if(grid1->opt_direction==GRID_3D_OPTIMIZE_RPT)
  {
    new_grid=create_3Dgrid_radial_major(npol, nrad, ntor);
  }
  else if(grid1->opt_direction==GRID_3D_OPTIMIZE_PRT)
  {
    new_grid=create_3Dgrid_poloidal_major(npol, nrad, ntor);
  }
  else
  {
    fprintf(stderr, "UNEXPECTED ERROR, UNKONWN Grid3DOptimization.\n");
    exit(EXIT_FAILURE);
  }
  #ifdef DEBUG
  printf("The size of new 3D grid: pol %5d rad %5d nt %5d.\n", npol, nrad, ntor);
  #endif
  //Merge the gird1 and grid2 in radial
  for(int k=0;k<ntor;k++)
  {
    for(int i=0;i<npol;i++)
    {
      for(int j=0;j<nrad1;j++)
      {
        set_point_3Dgrid(new_grid,i,j,k,
                         get_r_3Dgrid(grid1,i,j,k),
                         get_z_3Dgrid(grid1,i,j,k),
                         get_phi_3Dgrid(grid1,i,j,k));
      }
      for(int j=nrad1;j<nrad;j++)
      {
        int idx = j - nrad1 + overlap;
        set_point_3Dgrid(new_grid,i,j,k,
                         get_r_3Dgrid(grid2,i,idx,k),
                         get_z_3Dgrid(grid2,i,idx,k),
                         get_phi_3Dgrid(grid2,i,idx,k));
      }
    }
  }
  #ifdef DEBUG
  printf("Successfully merge two 3D grids along radial direction.\n");
  #endif
  return new_grid;
}


//reverse the coordinates in the poloidal direction.
void reverse_3Dgrid_poloidal_direction(ThreeDimGrid* grid);

//reverse the coordinates in the radial direction.
void reverse_3Dgrid_radial_direction(ThreeDimGrid* grid);

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

void write_EMC3_3Dgrid_to_EMC3_format(ThreeDimGrid* g, char* filename, bool reverse_pol, bool reverse_rad)
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
  const int np=g->npol;
  const int nr=g->nrad;
  const int nt=g->ntor;
  fprintf(fp, "%10d%10d%10d\n",nr, np, nt);
  for(int it=0; it<nt; it++)
  {
    double phi=get_phi_3Dgrid(g,0,0,it);
    fprintf(fp, "%15.8f\n",phi);

    for(int ip=0; ip<np; ip++)
    {
      const int ii = reverse_pol ? (np - 1 - ip) : ip;
      for(int ir=0; ir<nr; ir++)
      {
        const int jj = reverse_rad ? (nr - 1 - ir) : ir;
        double r=get_r_3Dgrid(g,ii,jj,it);
        r=100*r; //meter to centimeter
        fprintf(fp, "%.12e\n", r);
      }
    }
    for(int ip=0; ip<np; ip++)
    {
      const int ii = reverse_pol ? (np - 1 - ip) : ip;
      for(int ir=0; ir<nr; ir++)
      {
        const int jj = reverse_rad ? (nr - 1 - ir) : ir;
        double z=get_z_3Dgrid(g,ii,jj,it);
        z=100*z; //meter to centimeter
        fprintf(fp, "%.12e\n", z);
      }
    }
  }
  fclose(fp);
  #ifdef DEBUG
  printf("Successfully write EMC3 format 3D GRID: %s.\n",filename);
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
        r=r/100; //centimeter to meter
        set_r_3Dgrid(grid3d, ip,ir,it, r);
        set_phi_3Dgrid(grid3d, ip,ir,it, phi);
      }
    }
    for(int ip=0; ip<np; ip++)
    {
      for(int ir=0; ir<nr; ir++)
      {
        fscanf(fp, "%lf", &z);
        z=z/100; //centimeter to meter
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

void radial_mapping_check_test(ThreeDimGrid* g1, int idx_r1, int idx_p1_s, int idx_p1_e,
                               ThreeDimGrid* g2, int idx_r2, int idx_p2_s, int idx_p2_e)
{
  if(!g1 || !g2)
  {
    fprintf(stderr, "Error: NULL input to in RADIAL MAPPING CHECK.\n");
    exit(EXIT_FAILURE);
  }
  int np1=idx_p1_e-idx_p1_s+1;
  int np2=idx_p2_e-idx_p2_s+1;

  if(np1 != np2)
  {
    fprintf(stderr, "Error: The sizes in poloidal direction are not identical in RADIAL MAPPING CHECK.\n");
    exit(EXIT_FAILURE);
  }

  int nt1=g1->ntor;
  int nt2=g2->ntor;
  if(nt1!=nt2)
  {
    fprintf(stderr, "Error: The sizes in toroidal direction are not identical in RADIAL MAPPING CHECK.\n");
    exit(EXIT_FAILURE);
  }
  bool pass=true;
  for(int it=0;it<nt1;it++)
  {
    double phi1=get_phi_3Dgrid(g1,0,0,it);
    double phi2=get_phi_3Dgrid(g2,0,0,it);
    if(fabs(phi1-phi2)>EPSILON_10)
    {
      fprintf(stderr, "Error: the %d-th phi of g1 and g2 are not identical.\n",it);
      fprintf(stderr, "g1 phi: %6f, g2 phi: %6f.\n",phi1,phi2);
      exit(EXIT_FAILURE);
    }
    for(int ip=0;ip<np1;ip++)
    {
      double r1=get_r_3Dgrid(g1,idx_p1_s+ip,idx_r1,it);
      double z1=get_z_3Dgrid(g1,idx_p1_s+ip,idx_r1,it);

      double r2=get_r_3Dgrid(g2,idx_p2_s+ip,idx_r2,it);
      double z2=get_z_3Dgrid(g2,idx_p2_s+ip,idx_r2,it);

      if (fabs(r1 - r2) > EPSILON_10 || fabs(z1 - z2) > EPSILON_10) 
      {
        pass=false;
        printf("Phi=%6f, ip=%6d: g1 R=%18.12f Z=%18.12f | g2 R=%18.12f Z=%18.12f "
               "(dR=%g, dZ=%g)\n",
               phi1, ip, r1, z1, r2, z2, r1 - r2, z1 - z2);
      }
      else
      {
        //using g1 values to assign g2 values to ensure the exactly the same
        set_point_3Dgrid(g2,idx_p2_s+ip,idx_r2,it, r1,z1,phi1);
      }
    }
  }
  if(pass==false)
  {
    fprintf(stderr, "Error: Failed RADIAL MAPPING CHECK.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    printf("Successful pass the RADIAL MAPPING CHECK.\n");
  }
}

void poloidal_extend_for_toroidal_mapping_test(ThreeDimGrid* grid, int n_tracing_step,
                                              ode_function* func,ode_solver* solver)
{
  if(!grid || !func || !solver)
  {
    fprintf(stderr, "Error: NULL input to poloidal_extend_for_toroidal_mapping_test.\n");
    exit(EXIT_FAILURE);
  }

  if(n_tracing_step<1)
  {
    fprintf(stderr, "Error: As least 1 step for n_tracing_step.\n");
    exit(EXIT_FAILURE);
  }

  const int ntor = grid->ntor;  // toroidal dimension unchanged here
  const int nrad = grid->nrad;  // radial dimension unchanged here

  /******************
   *  Inner portion *
   ******************/
  // Extend one logical cell at poloidal head and tail (+2 total).
  // This only updates logical size/offset; it does NOT move data.
  expand_3Dgrid(grid, 1, 0, 0, 0, 0, 0);  // add head cell at pol=0
  expand_3Dgrid(grid, 0, 1, 0, 0, 0, 0);  // add tail cell at pol=grid->npol-1

  //restore magnetic field direction to ensure the correct
  restore_3D_mag_direction(func);
  reverse_3D_mag_direction(func);
  
  for(int ir=0;ir<nrad;ir++)
  {
    double pt[3];
    //BECAREFULL, 1 is the origanl start point before expansion,
    pt[0]=get_r_3Dgrid(grid,1,ir,ntor-1);
    pt[1]=get_z_3Dgrid(grid,1,ir,ntor-1);
    pt[2]=get_phi_3Dgrid(grid,1,ir,ntor-1);

    double pt_tmp[3] = { pt[0], pt[1], pt[2] };
    double next_pt_tmp[3];
    double t_tmp = 0.0;

    for(int ii=0;ii<n_tracing_step;ii++)
    {
      // Take one integration step
      solver->next_step(solver->step_size, &t_tmp, pt_tmp, next_pt_tmp, 
                        solver->solver_data, func);
      t_tmp=t_tmp+solver->step_size;

      pt_tmp[0] = next_pt_tmp[0];
      pt_tmp[1] = next_pt_tmp[1];
      pt_tmp[2] = next_pt_tmp[2];
    }

    //Check the point
    if(pt_tmp[2]<pt[2])
    {
      fprintf(stderr,"Unexpected results.\n");
      exit(EXIT_FAILURE);
    }

    #ifdef DEBUG
    printf("The %d-th radial at inner portrain expend to %.12f %.12f.\n", ir, pt_tmp[0],pt_tmp[1]);
    #endif


    for(int it=0;it<ntor;it++)
    {
      double phi=get_phi_3Dgrid(grid,1,ir,it);
      set_point_3Dgrid(grid,0,ir,it,pt_tmp[0],pt_tmp[1],phi);
    }
  }
  
  
  /******************
   * Outer portion  *
   ******************/
  
  //restore magnetic field direction to ensure the correct
  restore_3D_mag_direction(func);
  // After two expansions, the last original poloidal layer sits at index (npol - 2).
  const int idx_pol = grid->npol - 2;
  if (idx_pol < 1) {
    fprintf(stderr, "Error: idx_pol out of range after expansion (idx_pol=%d, npol=%d).\n",
            idx_pol, grid->npol);
    exit(EXIT_FAILURE);
  }
  for(int ir=0;ir<nrad;ir++)
  {
    double pt[3];
    //BECAREFUL, 1 is the origanl start point before expansion,
    
    pt[0]=get_r_3Dgrid(grid,idx_pol,ir,0);
    pt[1]=get_z_3Dgrid(grid,idx_pol,ir,0);
    pt[2]=get_phi_3Dgrid(grid,idx_pol,ir,0);

    double pt_tmp[3] = { pt[0], pt[1], pt[2] };
    double next_pt_tmp[3];
    double t_tmp = 0.0;

    for(int ii=0;ii<n_tracing_step;ii++)
    {
      // Take one integration step
      solver->next_step(solver->step_size, &t_tmp, pt_tmp, next_pt_tmp, 
                        solver->solver_data, func);
      t_tmp=t_tmp+solver->step_size;

      pt_tmp[0] = next_pt_tmp[0];
      pt_tmp[1] = next_pt_tmp[1];
      pt_tmp[2] = next_pt_tmp[2];
    }

    //Check the point
    if(pt_tmp[2]>pt[2])
    {
      fprintf(stderr,"Unexpected results.\n");
      exit(EXIT_FAILURE);
    }

    for(int it=0;it<ntor;it++)
    {
      double phi=get_phi_3Dgrid(grid,idx_pol,ir,it);
      set_point_3Dgrid(grid,idx_pol+1,ir,it,pt_tmp[0],pt_tmp[1],phi);
    }

    #ifdef DEBUG
    printf("The %d-th radial at inner portrain expend to %.12f %.12f.\n", ir, pt_tmp[0],pt_tmp[1]);
    #endif
  }
}
