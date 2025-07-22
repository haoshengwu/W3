#include "twodimgridgen.h"
#include "curve.h"
#include <math.h>
#include "mathbase.h"
#include "utils.h"
#include "datastructure.h"

#define _POSIX_C_SOURCE 200112L

#define DEFAULT_MARGIN 20
#define MAX_NUM_TRACING 40000

#ifndef EPS_TDGG
#define EPS_TDGG 1.0E-12
#endif

#ifndef NRELAX
#define NRELAX 5000
#endif

#ifndef RLCEPT
#define RLCEPT 1.0E-6
#endif

#ifndef RELAX
#define RELAX 0.2
#endif


/*
 * Merge an array of doubly-linked lists into a single Curve.
 * ──────────────────────────────────────────────────────────
 * Input
 *   list : array of DLListWithOptions, each containing a head node
 *          and a reverse flag that indicates traversal direction.
 *   n    : number of list segments in the array
 *
 * Rules
 *   1.  The first segment is copied in full.
 *   2.  For every subsequent segment, if its first node duplicates
 *       the last point that was added to the result, that first node
 *       is skipped.
 *   3.  The original linked lists are NOT modified.
 */
Curve* connect_DLList_for_curve(DLListWithOptions* list, int n)
{
  if (!list || n <= 0) {
    fprintf(stderr, "Error: invalid input to connect_DLList_for_curve.\n");
    exit(EXIT_FAILURE);
  }

  /* Allocate the result curve with an initial capacity. */
  Curve* result = create_curve(256);
  if (!result) {
    fprintf(stderr, "Error: failed to allocate result curve.\n");
    exit(EXIT_FAILURE);
  }

  for (int seg = 0; seg < n; ++seg) {

    /* Validate current segment head */
    DLListNode* head = list[seg].head;
    if (!head) {
      fprintf(stderr, "Error: DLList segment %d is NULL.\n", seg);
      exit(EXIT_FAILURE);
    }

    /* Choose traversal start node based on the reverse flag */
    DLListNode* node = list[seg].reverse
                       ? get_DLList_tailnode(head)   /* start from tail */
                       : head;                       /* start from head */

    /* Flag to indicate the very first node of this segment */
    int is_first_point = 1;

    /* Traverse the current linked list */
    while (node) {
      /* For segments after the first, skip the first node if it
         duplicates the last point already in the curve.          */
      if (seg > 0 && is_first_point && result->n_point > 0) {
        double last_x = result->points[result->n_point - 1].x;
        double last_y = result->points[result->n_point - 1].y;
        if (fabs(last_x - node->r) < 1e-10 &&
            fabs(last_y - node->z) < 1e-10) {
          /* Duplicate found -- do not add this point */
          node = list[seg].reverse ? node->prev : node->next;
          is_first_point = 0;
          continue;
        }
      }

      /* Add current node’s coordinates to the curve */
      add_last_point_curve(result, node->r, node->z);

      /* Advance */
      node = list[seg].reverse ? node->prev : node->next;
      is_first_point = 0;     /* We are past the first point now   */
    }
  }

  return result;
}


GridZone* create_sn_CARRE2D_GridZone(GridZoneInfo* gzinfo, SepDistStr* sepdist)
{
  // Allocate and copy name/topo
  GridZone* gz = malloc(sizeof(GridZone));
  if (!gz) 
  {
    fprintf(stderr, "Failed to allocate GridZone.\n");
    exit(EXIT_FAILURE);
  }

  int len = strlen(gzinfo->topo) + 1;
  gz->topo = malloc(len);
  strcpy(gz->topo, gzinfo->topo);

  len = strlen(gzinfo->name) + 1;
  gz->name = malloc(len);
  strcpy(gz->name, gzinfo->name);

  // Copy numerical arrays
  int nr = gzinfo->nr;
  gz->nr = nr;

  gz->start_point_R = malloc(nr * sizeof(double));
  gz->start_point_Z = malloc(nr * sizeof(double));
  gz->guard_start   = malloc(nr * sizeof(double));
  gz->guard_end     = malloc(nr * sizeof(double));
  gz->pasmin        = malloc(nr * sizeof(double));

  memcpy(gz->start_point_R, gzinfo->start_point_R, nr * sizeof(double));
  memcpy(gz->start_point_Z, gzinfo->start_point_Z, nr * sizeof(double));
  memcpy(gz->guard_start,   gzinfo->guard_start,   nr * sizeof(double));
  memcpy(gz->guard_end,     gzinfo->guard_end,     nr * sizeof(double));
  memcpy(gz->pasmin,        gzinfo->pasmin,        nr * sizeof(double));

  gz->end_curve=copy_curve(gzinfo->end_curve);


  // create the first boundary curve and grid point curve
  int n_segm1=gzinfo->n_polsegm1;
  CurveWithOptions*  option_gp_c= malloc(n_segm1*sizeof(CurveWithOptions));
  DLListWithOptions* option_c = malloc(n_segm1*sizeof(DLListWithOptions));
  for(int i=0; i<n_segm1; i++)
  {
    printf("DEBUG n_segm1 %d\n",n_segm1);
    int idx=sepdist->index[gzinfo->seplineidx1[i]];
    printf("DEBUG gzinfo->seplineidx1 %d\n",gzinfo->seplineidx1[i]);
    printf("DEBUG idx %d\n", idx);
    option_gp_c[i].curve=sepdist->edges[idx]->gridpoint_curve;
    option_gp_c[i].reverse=gzinfo->reverse_segm1[i];

    option_c[i].head=sepdist->edges[idx]->head;
    option_c[i].reverse=gzinfo->reverse_segm1[i];
   }
  gz->first_bnd=true;
  gz->first_gridpoint_curve=connect_curves_for_curve(option_gp_c,n_segm1);
  gz->first_bnd_curve=connect_DLList_for_curve(option_c, n_segm1);

  gz->sec_bnd=false;
  gz->sec_gridpoint_curve=NULL;
  gz->sec_bnd_curve=NULL;
  free(option_c);
  free(option_gp_c);
  return gz;
}

void free_GridZone(GridZone* gz)
{
  if (!gz) return;

  // Free strings
  if (gz->topo) {
    free(gz->topo);
    gz->topo = NULL;
  }

  if (gz->name) {
    free(gz->name);
    gz->name = NULL;
  }

  // Free double arrays
  if (gz->start_point_R) {
    free(gz->start_point_R);
    gz->start_point_R = NULL;
  }

  if (gz->start_point_Z) {
    free(gz->start_point_Z);
    gz->start_point_Z = NULL;
  }

  if (gz->guard_start) {
    free(gz->guard_start);
    gz->guard_start = NULL;
  }

  if (gz->guard_end) {
    free(gz->guard_end);
    gz->guard_end = NULL;
  }

  if (gz->pasmin) {
    free(gz->pasmin);
    gz->pasmin = NULL;
  }

  // Free curves if allocated
  if (gz->end_curve) {
    free_curve(gz->end_curve);
    gz->end_curve = NULL;
  }

  if (gz->first_bnd_curve) {
    free_curve(gz->first_bnd_curve);
    gz->first_bnd_curve = NULL;
  }

  if (gz->first_gridpoint_curve) {
    free_curve(gz->first_gridpoint_curve);
    gz->first_gridpoint_curve = NULL;
  }

  if (gz->sec_bnd_curve) {
    free_curve(gz->sec_bnd_curve);
    gz->sec_bnd_curve = NULL;
  }

  if (gz->sec_gridpoint_curve) {
    free_curve(gz->sec_gridpoint_curve);
    gz->sec_gridpoint_curve = NULL;
  }
  // Free the GridZone structure itself
  free(gz);
}

typedef struct
{
/*********************************
* Input parameters
**********************************/
    int np; //poloidal point number
  
    //previous curve which is coresponding to (ir-1) gridcurve. 
    Curve *prev_c;

    //curve points which is coresponding to (ir )gridcurve.
    Curve *curr_c;
    
    //The length distributions of the grid pionts in the previous gridcurve. 
    double *len_prev_gpt_c;

    //the points' coordiantes for the grid which are along the prev_curve.
    Curve *prev_gpt_c;

    //guardlength at the start for curve
    double guard_top;
    //guardlength at the end for curve
    double guard_end;
    //minimum distance between mesh points in the curve
    double pasmin;
/*********************************
* Output parameters
**********************************/
    //The length of the mesh pionts in the curve, is the length from the start point to the point
    double *len_curr_gpt_c;
    Curve *curr_gpt_c;

} GirdTubeStr;

//Only gether the entries together, not create new one!!! and NO FREE function.
static GirdTubeStr* create_GridTube(Curve *prev_c,
                                    Curve *prev_gpt_c,
                                    double *len_prev_gpt_c, 
                                    Curve *curr_c, 
                                    Curve *curr_gpt_c, //output
                                    double *len_curr_gpt_c,//output
                                    double guard_top, 
                                    double guard_end, 
                                    double pasmin
                                    )
{
  if(guard_end<0.0||guard_end<0.0)
  {
    printf("WARNING: Guard length is less than ZERO!\n");
    printf("WARNING: Guard length is set to 0.0!\n");
    guard_end=0.0;
    guard_top=0.0;
  }
  if(!prev_c||!prev_gpt_c||!len_prev_gpt_c||!curr_c||!curr_gpt_c||!len_curr_gpt_c)
  {
    fprintf(stderr,"Empty input for create_GridTube.\n");
    exit(EXIT_FAILURE);
  }
  if(prev_gpt_c->n_point!=curr_gpt_c->n_point)
  {
    fprintf(stderr,"The size of prev_gpt_c is not consistent with curr_gpt_c.\n");
    exit(EXIT_FAILURE);
  }
  GirdTubeStr* gridtube = malloc(sizeof(GirdTubeStr));
  if(!gridtube)
  {
    fprintf(stderr,"Failed to allocate gridtube.\n");
    exit(EXIT_FAILURE);
  }
  gridtube->np=prev_gpt_c->n_point;
  gridtube->pasmin=pasmin;
  gridtube->guard_top=guard_top;
  gridtube->guard_end=guard_end;
  gridtube->prev_c=prev_c;
  gridtube->curr_c=curr_c;
  gridtube->prev_gpt_c=prev_gpt_c;
  gridtube->curr_gpt_c=curr_gpt_c;
  gridtube->len_prev_gpt_c= len_prev_gpt_c;
  gridtube->len_curr_gpt_c= len_curr_gpt_c;
  return gridtube;
}

/*====================================================================
| Adaptive 2D Grid System with Storage Order Optimization
| Just choose the optimization direction when creating the grid;
| all implementation details are fully encapsulated.
|
| This version uses 32-byte aligned memory allocation for GridPoint
| arrays, suitable for SIMD/high-performance scenarios.
|====================================================================*/

TwoDimGrid* create_2Dgrid_poloidal_major(int npol, int nrad)
{
    return create_2Dgrid_optimized_for(npol, nrad, GRID_OPTIMIZE_FOR_IP);
}

TwoDimGrid* create_2Dgrid_radial_major(int npol, int nrad)
{
    return create_2Dgrid_optimized_for(npol, nrad, GRID_OPTIMIZE_FOR_IR);
}

TwoDimGrid* create_2Dgrid_optimized_for(int npol, int nrad, GridOptimization opt)
{
    int margin_pol = DEFAULT_MARGIN;
    int margin_rad = DEFAULT_MARGIN;

    TwoDimGrid* grid = malloc(sizeof(TwoDimGrid));
    if (!grid) {
        fprintf(stderr, "Failed to allocate TwoDimGrid.\n");
        exit(EXIT_FAILURE);
    }

    grid->npol = npol;
    grid->nrad = nrad;
    grid->cap_npol = npol + 2 * margin_pol;
    grid->cap_nrad = nrad + 2 * margin_rad;
    grid->offset_pol = margin_pol;
    grid->offset_rad = margin_rad;
    grid->opt_direction = opt;

    // Allocate aligned memory for high performance (32 bytes for AVX/SIMD)
    size_t total = (size_t)grid->cap_npol * grid->cap_nrad;
    size_t size = total * sizeof(GridPoint);
    
    // use aligned_alloc (C11)
    grid->points = aligned_alloc(32, size);

    if (!grid->points) {
        fprintf(stderr, "Failed to allocate aligned grid points.\n");
        free(grid);
        exit(EXIT_FAILURE);
    }
    // Zero-initialize all points (calloc-style effect)
    for (size_t i = 0; i < total; ++i) {
        grid->points[i].x = 0.0;
        grid->points[i].y = 0.0;
    }
    return grid;
}

void free_2Dgrid(TwoDimGrid* grid)
{
    if (!grid) return;
    free(grid->points);
    grid->points = NULL;
    free(grid);
}

static inline size_t idx_2Dgrid(const TwoDimGrid* g, int ip, int ir)
{
    if (g->opt_direction == GRID_OPTIMIZE_FOR_IP) {
        // Row-major: optimize ip access
        return (size_t)(g->offset_rad + ir) * g->cap_npol
             + (size_t)(g->offset_pol + ip);
    } else {
        // Column-major: optimize ir access
        return (size_t)(g->offset_pol + ip) * g->cap_nrad
             + (size_t)(g->offset_rad + ir);
    }
}

static inline bool is_valid_index(const TwoDimGrid* g, int ip, int ir) 
{
    return ip >= 0 && ip < g->npol && ir >= 0 && ir < g->nrad;
}

GridPoint* get_point_2Dgrid(const TwoDimGrid* g, int ip, int ir)
{
    return &g->points[idx_2Dgrid(g, ip, ir)];
}
double get_x_2Dgrid(const TwoDimGrid* g, int ip, int ir)
{
    return g->points[idx_2Dgrid(g, ip, ir)].x;
}

double get_y_2Dgrid(const TwoDimGrid* g, int ip, int ir)
{
    return g->points[idx_2Dgrid(g, ip, ir)].y;
}

void set_point_2Dgrid(TwoDimGrid* g, int ip, int ir, double x, double y)
{
    GridPoint* p = &g->points[idx_2Dgrid(g, ip, ir)];
    p->x = x;
    p->y = y;
}

int get_npol_2Dgrid(const TwoDimGrid* g) { return g->npol; }
int get_nrad_2Dgrid(const TwoDimGrid* g) { return g->nrad; }

/*--------------------------------------------------------------------
  High-performance traversal - automatically uses optimal order
--------------------------------------------------------------------*/
void foreach_point_2Dgrid(const TwoDimGrid* g, 
                          void (*callback)(int ip, int ir, double x, double y, void* userdata),
                          void* userdata)
{
    if (g->opt_direction == GRID_OPTIMIZE_FOR_IP) {
        // Optimize ip: outer loop ir, inner loop ip
        for (int ir = 0; ir < g->nrad; ir++) {
            for (int ip = 0; ip < g->npol; ip++) {
                double x = get_x_2Dgrid(g, ip, ir);
                double y = get_y_2Dgrid(g, ip, ir);
                callback(ip, ir, x, y, userdata);
            }
        }
    } else {
        // Optimize ir: outer loop ip, inner loop ir
        for (int ip = 0; ip < g->npol; ip++) {
            for (int ir = 0; ir < g->nrad; ir++) {
                double x = get_x_2Dgrid(g, ip, ir);
                double y = get_y_2Dgrid(g, ip, ir);
                callback(ip, ir, x, y, userdata);
            }
        }
    }
}

void expand_2Dgrid(TwoDimGrid* g, 
                   int add_pol_head, int add_pol_tail,
                   int add_rad_head, int add_rad_tail)
{
  if (!g) return;
  if (add_pol_head < 0 || add_pol_tail < 0 || add_rad_head < 0 || add_rad_tail < 0)
  {
    fprintf(stderr,"expaned number should>=0 for expand_2Dgrid");
    exit(EXIT_FAILURE);
  }
  
  int new_npol = g->npol + add_pol_head + add_pol_tail;
  int new_nrad = g->nrad + add_rad_head + add_rad_tail;
  int new_offset_pol = g->offset_pol - add_pol_head;
  int new_offset_rad = g->offset_rad - add_rad_head;
  
  if (new_offset_pol < 0 || new_offset_rad < 0 ||
      (new_offset_pol + new_npol) > g->cap_npol ||
      (new_offset_rad + new_nrad) > g->cap_nrad) 
  {
    fprintf(stderr,
           "expand_2Dgrid: expansion exceeds allocated memory (capacity: pol=%d, rad=%d). "
           "Requested offset_pol=%d, npol=%d, offset_rad=%d, nrad=%d\n",
            g->cap_npol, g->cap_nrad,
            new_offset_pol, new_npol, new_offset_rad, new_nrad);
    fprintf(stderr,"In the future, dynamic expansion will be supported.\n");
    exit(EXIT_FAILURE); 
  }
  g->offset_pol = new_offset_pol;
  g->offset_rad = new_offset_rad;
  g->npol = new_npol;
  g->nrad = new_nrad;
}

static void write_2Dgrid_callback(int ip, int ir, double x, double y, void* userdata) {
    FILE* f = (FILE*)userdata;
    fprintf(f, "%.12f %.12f\n", x, y);
}

void write_2Dgrid(const TwoDimGrid* g, char* filename)
{
  if (!filename || !g) 
  {
    fprintf(stderr, "Error: NULL input to write_2Dgrid.\n");
    exit(EXIT_FAILURE);
  }

  FILE *fp = fopen(filename, "w");
  if (!fp) 
  {
    fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "# %d %d %s\n",g->npol, g->nrad,
         (g->opt_direction == GRID_OPTIMIZE_FOR_IP) ? "GRID_OPTIMIZE_FOR_IP" : "GRID_OPTIMIZE_FOR_IR");
  foreach_point_2Dgrid(g, write_2Dgrid_callback, fp);
  printf("Finish writing 2Dgrid to %s file.\n",filename);
  fclose(fp);
}


TwoDimGrid* load_2Dgrid_from_file(char* filename)
{
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    perror("Failed to open input file");
    exit(EXIT_FAILURE);
  }
  int npol, nrad;
  char opt_str[32];
  if (fscanf(fp, "# %d %d %s\n", &npol, &nrad, opt_str) != 3)
  {
    fprintf(stderr, "Error: invalid header format in file \"%s\"\n", filename);
    fclose(fp);
    exit(EXIT_FAILURE);
  }
  GridOptimization opt_direction;
  TwoDimGrid* grid;
  if (strcmp(opt_str, "GRID_OPTIMIZE_FOR_IP") == 0)
  {
    opt_direction = GRID_OPTIMIZE_FOR_IP;
    grid=create_2Dgrid_poloidal_major(npol,nrad);
  }
  else if (strcmp(opt_str, "GRID_OPTIMIZE_FOR_IR") == 0)
  {
    opt_direction = GRID_OPTIMIZE_FOR_IR;
    grid=create_2Dgrid_radial_major(npol,nrad);
  }
  else
  {
    fprintf(stderr, "Error: unknown optimization direction \"%s\"\n", opt_str);
    fclose(fp);
    exit(EXIT_FAILURE);
  }
  // Use the same traversal order based on optimization direction
  if (opt_direction == GRID_OPTIMIZE_FOR_IP)
  {
    // Same as write: outer loop ir, inner loop ip
    for (int ir = 0; ir < nrad; ir++) 
    {
      for (int ip = 0; ip < npol; ip++) 
      {
        double x, y;
        if (fscanf(fp, "%lf %lf", &x, &y) != 2)
        {
          fprintf(stderr, "Error: failed to read grid point (%d, %d)\n", ip, ir);
          fclose(fp);
          free_2Dgrid(grid);
          exit(EXIT_FAILURE);
        }
        set_point_2Dgrid(grid, ip, ir, x, y);
      }
    }
  } 
  else if (opt_direction == GRID_OPTIMIZE_FOR_IR)
  {
    // Same as write: outer loop ip, inner loop ir
    for (int ip = 0; ip < npol; ip++) 
    {
      for (int ir = 0; ir < nrad; ir++) 
      {
        double x, y;
        if (fscanf(fp, "%lf %lf", &x, &y) != 2)
        {
          fprintf(stderr, "Error: failed to read grid point (%d, %d)\n", ip, ir);
          fclose(fp);
          free_2Dgrid(grid);
          exit(EXIT_FAILURE);
        }
        set_point_2Dgrid(grid, ip, ir, x, y);
      }
    }
  }
  else
  {
    fprintf(stderr, "Error: unknown optimization direction \"%s\"\n", opt_str);
    fclose(fp);
    exit(EXIT_FAILURE);
  }

  printf("Finish loading 2Dgrid from %s file.\n", filename);
  fclose(fp);
  
  return grid;
}



/****************************************************************************
*    INHERIT FROM CAREE
****************************************************************************/

static double calc_length(CurvePoint* p1, CurvePoint* p2) 
{
    return hypot(p1->x-p2->x, p1->y-p2->y);
}


static double cosine_term(CurvePoint* base, CurvePoint* p1, CurvePoint* p2) 
{
    double x1 = p1->x-base->x;
    double x2 = p2->x-base->x;
    double y1 = p1->y-base->y;
    double y2 = p2->y-base->y;
    return dot_product(x1,y1,x2,y2);
}

typedef struct
{
  size_t n;
  double *ort;
  double *ortpur;
  double *propo;
  double *varr;
  double *tot;
} CarreOrthoStr;

static CarreOrthoStr* allocate_CarreOrthoProp(size_t n)
{
  if(n<2)
  {
    fprintf(stderr,"Empty input for CarreOrthoStr.\n");
    exit(EXIT_FAILURE);
  }
  CarreOrthoStr* orthogonal=malloc(sizeof(CarreOrthoStr));
  if(!orthogonal)
  {
    fprintf(stderr,"Failed to allocate CarreOrthoStr.\n");
    exit(EXIT_FAILURE);
  }
  orthogonal->n = n;
  orthogonal->ort = (double *)calloc(n, sizeof(double));
  orthogonal->ortpur = (double *)calloc(n, sizeof(double));
  orthogonal->propo  = (double *)calloc(n, sizeof(double));
  orthogonal->varr = (double *)calloc(n, sizeof(double));
  orthogonal->tot = (double *)calloc(n, sizeof(double));
  return orthogonal;
}

static void free_CarreOrthoProp(CarreOrthoStr *orthogonal)
{
  if(!orthogonal) return;
  free(orthogonal->ort);orthogonal->ort=NULL;
  free(orthogonal->ortpur); orthogonal->ortpur=NULL;
  free(orthogonal->propo); orthogonal->propo=NULL;
  free(orthogonal->varr); orthogonal->varr=NULL;
  free(orthogonal->tot); orthogonal->tot=NULL;
  free(orthogonal);
}

/*******************************************************************************
* This function calculate the orthogonalirty of the mesh points on two mesh curves.
* This fuction is refered to Fortran code CARRE/clort.F .
********************************************************************************/
//len_prev_gpt_c is length distribution for previous grid point curve, 0.0 for the 1st and total length for the last
//Curve *prev_gpt_c is the previous curve contain the grid point
//len_curr_gpt_c is length distribution for current  grid point curve, 0.0 for the 1st and total length for the last
//Curve *curr_gpt_c is the current curve contain the grid point
//'current curve' means we want the grid point in the current cure based on previous curve.

static void calc_ortho_from_CARRE(int n_point,
                                  const double* len_prev_gpt_c,
                                  const Curve* prev_gpt_c,
                                  double* len_curr_gpt_c,
                                  const Curve* curr_gpt_c,
                                  double guard_top,
                                  double guard_end,
                                  double pasmin,
                                  CarreOrthoStr* orthogonal)

{
  if (!prev_gpt_c || !curr_gpt_c || !orthogonal) 
  {
    fprintf(stderr, "Null input in calc_ortho_CARRE_explicit.\n");
    exit(EXIT_FAILURE);
  }

  if (prev_gpt_c->n_point != n_point || curr_gpt_c->n_point != n_point) 
  {
    fprintf(stderr, "The number of points is inconsistent in the curves!\n");
    exit(EXIT_FAILURE);
  }

/********************************************
 * 
 * i-1 \        i+1/ previous gridcurve
 *      \         /
 *       ----i----
 *           |
 *    -------I--------
 *   /                \
 *  / I-1           I+1\  current gridcurve
 *
 *********************************************/
  const double un = 1.0;

  double l1p, l2p, l12t;
  l1p  = calc_length(&prev_gpt_c->points[0],
                     &prev_gpt_c->points[1]);
  l2p  = calc_length(&curr_gpt_c->points[0],
                     &curr_gpt_c->points[1]);
  l12t = calc_length(&prev_gpt_c->points[0],
                     &curr_gpt_c->points[0]);

  double g1 = guard_top;
  double g2 = guard_end;

  for (int i = 1; i < n_point - 1; i++) 
  {
    double l12;
    if (g1 > EPS_TDGG && g2 > EPS_TDGG) 
    {
      double fac1 = pow(g1 / (g1 + len_prev_gpt_c[i]), 2.0);
      double fac2 = pow(g2 / (g2 + len_prev_gpt_c[n_point - 1] -
                              len_prev_gpt_c[i]), 2.0);
      double fac  = fac1 * (un - fac2) + fac2;
      l12 = fac * len_prev_gpt_c[n_point - 1] +
            (1.0 - fac) * l12t;
    } 
    // Modified for different situations.
    else if (g1 > EPS_TDGG && g2 < EPS_TDGG)
    {
      double fac1 =  pow(g1 / (g1 + len_prev_gpt_c[i]), 2.0);
      double fac = fac1;
      l12 = fac * len_prev_gpt_c[n_point - 1] +
            (1.0 - fac) * l12t;
    }
    else if (g1 < EPS_TDGG && g2 > EPS_TDGG)
    {
      double fac2 = pow(g2 / (g2 + len_prev_gpt_c[n_point - 1] -
                  len_prev_gpt_c[i]), 2.0);
      double fac = fac2;
      l12 = fac * len_prev_gpt_c[n_point - 1] +
            (1.0 - fac) * l12t;
    }
    else
    {
      l12 = l12t;
    }

    double l1m = l1p;
    double l2m = l2p;

    l1p = calc_length(&prev_gpt_c->points[i],
                      &prev_gpt_c->points[i + 1]);
    l2p = calc_length(&curr_gpt_c->points[i],
                      &curr_gpt_c->points[i + 1]);

    double cs1 = 0.0, cs2 = 0.0, cs3 = 0.0, cs4 = 0.0;
    if (i - 1 >= 0 && i + 1 < n_point) {
      cs1 = cosine_term(&curr_gpt_c->points[i],
                        &prev_gpt_c->points[i],
                        &curr_gpt_c->points[i - 1]);
      cs2 = cosine_term(&curr_gpt_c->points[i],
                        &prev_gpt_c->points[i],
                        &curr_gpt_c->points[i + 1]);
      cs3 = cosine_term(&prev_gpt_c->points[i],
                        &curr_gpt_c->points[i],
                        &curr_gpt_c->points[i - 1]);
      cs4 = cosine_term(&prev_gpt_c->points[i],
                        &curr_gpt_c->points[i],
                        &curr_gpt_c->points[i + 1]);
    } else 
    {
      printf("Warning: in calc_ortho_CARRE, Out of bounds access at i=%d\n", i);
    }

    cs1 /= (l2m * l12);
    cs2 /= (l2p * l12);
    cs3 /= (l1m * l12);
    cs4 /= (l1p * l12);

    double f1 = cs2 + cs3 - cs1 - cs4;
    double f2 = -(pow(g1 / len_curr_gpt_c[i], 2.0) +
                  pow(g2 / (len_curr_gpt_c[n_point - 1] -
                            len_curr_gpt_c[i]), 2.0)) *
                (len_curr_gpt_c[i] -
                 len_prev_gpt_c[i] /
                 len_prev_gpt_c[n_point - 1] *
                 len_curr_gpt_c[n_point - 1]) /
                (pasmin + g1 + g2);
    double f3 = pow(pasmin /
                    (len_curr_gpt_c[i] -
                     len_curr_gpt_c[i - 1]), 2.0) -
                pow(pasmin /
                    (len_curr_gpt_c[i + 1] -
                     len_curr_gpt_c[i]), 2.0);
    // printf("DEBUG f1 %lf f2 %lf f3 %lf\n", f1,f2,f3);
    orthogonal->ort[i]    = f1 + f2 + f3;
    orthogonal->ortpur[i] = f1;
    orthogonal->propo[i]  = f2;
    orthogonal->varr[i]   = f3;
    orthogonal->tot[i]    = f1 + f2 + f3;
  }
}

/*******************************************************************************
* This function calculate the mesh points in the curve which have a good orthogonalirty
* using a secant method according to R. Marchand Computer Physics Communications, 1996, 96.2-3: 232-246.
* In carre, the guard length is depend on separatrix. In this code, it dependt on previous cureve.
* Then the guard length are vertorized(different curves have different value). 
********************************************************************************/

void calc_points_from_CARRE(GirdTubeStr *tube)
{
  size_t n_point = tube->np;
  CarreOrthoStr* tmp_ortho= allocate_CarreOrthoProp(n_point);
  CarreOrthoStr* ortho= allocate_CarreOrthoProp(n_point);
  
  if(NRELAX >= 0)
  {
    double d1 = 0.0;
    size_t ipol1 = 1;
    size_t ipoln = n_point - 1;
//1. we first arrange the points proportionally to those of the previous line.
    //tmp_length_points is coresponding to l1 in mailrg.F
    //tube->length_points is coresponding to l2 mailrg.F
    double *tmp_length_points = (double *)calloc(n_point, sizeof(double));

    //initialize all points to zero;
    Curve* tmp_gpt_c=create_curve(n_point);

    expand_curve_size_with_NaN(tmp_gpt_c,n_point);
    tmp_length_points[0] = 0.0; 
    double length = total_length_curve(tube->curr_c);
    tmp_length_points[n_point-1] = length;

    double prev_length = tube->len_prev_gpt_c[n_point-1];
    
    // the first and the last point of mesh points are the same with curve
    int idx = tube->curr_c->n_point-1;
    tube->curr_gpt_c->points[0].x=tube->curr_c->points[0].x;
    tube->curr_gpt_c->points[0].y=tube->curr_c->points[0].y;
    tube->curr_gpt_c->points[n_point-1].x=tube->curr_c->points[idx].x;
    tube->curr_gpt_c->points[n_point-1].y=tube->curr_c->points[idx].y;

    tmp_gpt_c->points[0].x = tube->curr_c->points[0].x;
    tmp_gpt_c->points[0].y = tube->curr_c->points[0].y;
    tmp_gpt_c->points[n_point-1].x = tube->curr_c->points[idx].x;
    tmp_gpt_c->points[n_point-1].y = tube->curr_c->points[idx].y;

    //DEBUG
    // write_curve("DEBUG_prev_gpt_c", tube->prev_gpt_c);
    // write_curve("DEBUG_prev_c",tube->prev_c);
    //printf("DEBUG length: %lf\n",length);
//    printf("debug in calc_points_CARRE line 149\n");

    // the remaing point from 1 to n_point-2
    for(int ipol = ipol1; ipol<ipoln; ipol++)
    {
      d1 = ruban_curve(tube->prev_c, &(tube->prev_gpt_c->points[ipol]),d1);

      // printf("DEBUG len_prev_gpt_c/prev_length %lf\n",tube->len_prev_gpt_c[ipol]/prev_length);

      tmp_length_points[ipol] = (tube->len_prev_gpt_c[ipol]/prev_length)*length;
      coordnates_in_curve(tube->curr_c, tmp_length_points[ipol], &(tmp_gpt_c->points[ipol]));
    }
    //debug******************************************
    // write_array(tube->len_prev_gpt_c, n_point, "DEBUG_len_prev_gpt_c");
    // write_array(tmp_length_points, n_point, "DEBUG_tmp_length_points");

    // write_curve("DEBUG_tube_curr_c",tube->curr_c);
    // write_curve("DEBUG_tmp_gpt_c", tmp_gpt_c);
    //***********************************************************************

 //   printf("debug in calc_points_CARRE line 159\n");

    if(NRELAX > 0)
    {
//2. we initialize the function which must be zero for an orthogonal distribution.
      calc_ortho_from_CARRE(n_point, 
                            tube->len_prev_gpt_c, tube->prev_gpt_c,
                            tmp_length_points, tmp_gpt_c,
                            tube->guard_top, tube->guard_end, tube->pasmin,
                            tmp_ortho);
      // write_array(tmp_ortho->ort,n_point, "DEBUG_tmp_ort");

//3. we proceed to a first displacement of the nodes
      tube->len_curr_gpt_c[0] = 0;
      tube->len_curr_gpt_c[n_point-1] = length;
  
      for(int ipol = 1; ipol<ipoln; ipol++)
      {
        if(tmp_ortho->ort[ipol] > 0.0)
        {
          tube->len_curr_gpt_c[ipol] = 0.9*tmp_length_points[ipol] + 0.1*tmp_length_points[ipol+1];
        }
        else
        {
          tube->len_curr_gpt_c[ipol] = 0.9*tmp_length_points[ipol] + 0.1*tmp_length_points[ipol-1];
        }
        // coord_CARRE(tube->curve,tube->n_curve,tube->length_points[ipol], tube->point_coord[ipol]);

        coordnates_in_curve(tube->curr_c, tube->len_curr_gpt_c[ipol], &(tube->curr_gpt_c->points[ipol]));
      }

    //debug******************************************
    // write_curve("first_replace", tube->curr_gpt_c);
    //************************************************************
    //todo: store the ortho value for the whole radregion
    // somort(ir)= somort(ir)+ (ort1(ipol)/nppol)
    // somortpur(ir)= somortpur (ir)+ (ortpur1(ipol)/nppol)
    // sompropo(ir)= sompropo(ir)+ (propo1(ipol)/nppol)
    // somvarr(ir)= somvarr(ir)+ (varr1(ipol)/nppol)
    // somtot(ir)= somtot(ir)+ (tot1(ipol)/nppol)
    
//4. we relax the points iteratively to achieve the best possible orthogonality
      double ortmax = 0.0;
      
      for (int i = 0; i<NRELAX; i++)
      {
        calc_ortho_from_CARRE(n_point, 
                            tube->len_prev_gpt_c, tube->prev_gpt_c,
                            tube->len_curr_gpt_c, tube->curr_gpt_c,
                            tube->guard_top, tube->guard_end, tube->pasmin,
                            ortho);
        ortmax = 0.0;

        for (int ipol=ipol1; ipol<ipoln; ipol++)
        {
          double ortho_current = ortho->ort[ipol];
          double ortho_diff = ortho_current - tmp_ortho->ort[ipol];
          double length_diff = tube->len_curr_gpt_c[ipol] - tmp_length_points[ipol];
          double del = 0.0;
          double ortho_diff_abs = fabs(ortho_diff);

          if (fabs(ortho_current) <= RLCEPT)
          {
            continue;
          }

          if(ortho_diff_abs > RLCEPT*RLCEPT)
          {
            del = -ortho_current*length_diff/ortho_diff;
          }
          else
          {
            del = 0.0;
          } 

          if(del > 0.0) 
          {
            del = min(del, RELAX*(tube->len_curr_gpt_c[ipol+1]-tube->len_curr_gpt_c[ipol]));
          }
          else
          {
            del = max(del, RELAX*(tube->len_curr_gpt_c[ipol-1]-tube->len_curr_gpt_c[ipol]));
          }

          if(fabs(del) > RLCEPT*RLCEPT)
          {
            tmp_length_points[ipol] = tube->len_curr_gpt_c[ipol];
            tmp_ortho->ort[ipol] = ortho->ort[ipol];
            tube->len_curr_gpt_c[ipol] = tmp_length_points[ipol]+del;
          }
          // printf("debug in calc_points_CARRE line 241\n"); 
          //printf("del: %.12f\n", del);
          //printf("ipol: %d, length_point %.12f\n",ipol, tube->length_points[ipol]);
          coordnates_in_curve(tube->curr_c, tube->len_curr_gpt_c[ipol], &(tube->curr_gpt_c->points[ipol]));

          ortmax=max(ortmax,fabs(ortho_current));
        }
        // printf("debug: ortmax %.12f\n",ortmax);
        // printf("i: %d\n",i);
        if ((ortmax < RLCEPT))
        {
          printf("Finish optimized after %d iterations.\n",i);
          break;
        }
      }
      if (ortmax > RLCEPT)
      {
        printf("WARING: After optimization, the orthognoality is not good\n");
        printf("please adjust the parameters\n");
      }
      #ifdef DEBUG
        printf("debug: after optimaztion ortmax %.12f\n",ortmax);
      #endif
      //TODO write the orthognonal values
      //   do ipol=ipol1,ipoln
      //   somortp(ir)= somortp(ir)+ (ort2(ipol)/nppol)
      //   somortpurp(ir)= somortpurp(ir)+ (ortpur2(ipol)/nppol)
      //   sompropop(ir)= sompropop(ir)+ (propo2(ipol)/nppol)
      //   somvarrp(ir)= somvarrp(ir)+ (varr2(ipol)/nppol)
      //   somtotp(ir)= somtotp(ir)+ (tot2(ipol)/nppol)
      //   enddo
      free(tmp_length_points);
      free_curve(tmp_gpt_c);
    }
    else
    {
      //ToDo: need to do similar things as CARRE.
      printf("not support NRELAX=0\n");
    }
    // write_curve("DEBUG_curr_gpt_c", tube->curr_gpt_c);

  }
  free_CarreOrthoProp(tmp_ortho);
  free_CarreOrthoProp(ortho);
}


// Used to determind the tracing direction. Because the magnetic field line may not consistent with our assumed diretion.
//Compare the initial direction of the magnetic field line with the direction of the poloidal boundary curve 
//using a cross product, and then determine whether they align with the expected direction specified for SOL, PFR, or CORE.
//ONLY USED FOR CARRE 2D grid
static void check_poloidal_direction(GridZone* gridzone, ode_function* func, ode_solver* solver)
{
  //Recover to initial 1.0;
  if(func->ndim==2)
  {
    func->rescale[0]=1.0;
    func->rescale[1]=1.0;
  }
  else
  {
    fprintf(stderr,"Unexpect error, the dimention should be 2.\n");
    exit(EXIT_FAILURE);
  }
  double dir;
  int start;
  if(strncmp(gridzone->name,"SOL",3)==0) 
  {
    dir=1.0;
    start=0;
  }
  else if (strncmp(gridzone->name,"PFR",3)==0)
  {
    dir=-1.0;
    start=0;
  }
  else if (strncmp(gridzone->name,"CORE",4)==0)
  {
    dir=-1.0;
    start=1;
  }
  else
  {
    fprintf(stderr,"Unexpect error, Please check the type of zone.\n");
    exit(EXIT_FAILURE);
  }

  if (gridzone->nr < 2 || gridzone->first_bnd_curve->n_point < start + 2) 
  {
  fprintf(stderr, "Insufficient points in GridZone to determine direction.\n");
  exit(EXIT_FAILURE);
  }

  double p1[2];
  p1[0]=gridzone->start_point_R[1];
  p1[1]=gridzone->start_point_Z[1];

  double p2[2];
  double t=0.00;
  solver->next_step(solver->step_size, &t, p1, p2, solver->solver_data, func);

  //???
  double p3[2]={gridzone->start_point_R[start], gridzone->start_point_Z[start]};
  double p4[2]={gridzone->start_point_R[start+1], gridzone->start_point_Z[start+1]};

  if(cross_product(p2[0]-p1[0], p2[1]-p1[1], p4[0]-p3[0], p4[1]-p3[1])*dir<0.0)
  {
    func->rescale[0]=-1.0;
    func->rescale[1]=-1.0;
    printf("DEBUG reverse the poloidal diretction in %s.\n",gridzone->name);
  }
  else
  {
    printf("DEBUG the poloidal diretction is correct in %s.\n",gridzone->name);
  }
}

// tracing from a grid zone start point and finally arrive to the end curve.
// Here we assume the tracing direction HAS BEEN CORRECTED by change_poloidal_direction
// We also asumme that the tracing line will arrive at the end curve.
// after use, do not forget to FREE it!!!!!
// We also asumme there is only ONE intersection on the end curve.
// THIS IS EXPECTED FOR 2D MAGNETIC FIELD. 
static Curve* create_2d_GridTubeCurve_2d_tracing(double start_R, double start_Z, Curve* end_curve, 
                                              ode_function* func,ode_solver* solver)
{
  if(func->ndim!=2)
  {
    fprintf(stderr, "create_2d_GridTubeCurve_2d_tracing used for 2D magnetic field.\n");
    exit(EXIT_FAILURE);
  }

  if (end_curve->n_point < 2) 
  {
    fprintf(stderr, "Error: end_curve has too few points.\n");
    exit(EXIT_FAILURE);
  }
  Curve* new_c=create_curve(MAX_NUM_TRACING);
  add_last_point_curve(new_c, start_R, start_Z);
  double curr_p[2]={start_R, start_Z};
  double next_p[2];
  double t=0.00;
  double step_size=0.1;
  while(true)
  {
    solver->next_step(step_size, &t, curr_p, next_p, solver->solver_data, func);
    bool intersection=false;
    bool isclosed=false;
    for(int i=1; i<end_curve->n_point;i++)
    {
      // return 0 means found intersection
      // new_c->n_point>2 otherwise for the core region the first element of curve is always
      // has intersection with end curve.
      if(has_intersection(curr_p[0],curr_p[1],next_p[0],next_p[1],
                          end_curve->points[i-1].x, end_curve->points[i-1].y,
                          end_curve->points[i].x, end_curve->points[i].y)==0
                          && new_c->n_point>2) 
      {
        double intsect_x, intsect_y;
        get_intersection_point(curr_p[0],curr_p[1],next_p[0],next_p[1],
                               end_curve->points[i-1].x, end_curve->points[i-1].y,
                               end_curve->points[i].x, end_curve->points[i].y,
                               &intsect_x, &intsect_y);
        if(fabs(intsect_x-start_R)<1.0E-8&&fabs(intsect_y-start_Z)<1.0E-8)
        {
          add_last_point_curve(new_c,start_R, start_Z);
        }
        else
        {
          add_last_point_curve(new_c, intsect_x, intsect_y);
        }
        intersection=true;
        break;
      }
    }
    if(intersection||isclosed)
    {
      break;
    }
    else
    {
      add_last_point_curve(new_c, next_p[0], next_p[1]);
      curr_p[0]=next_p[0];
      curr_p[1]=next_p[1];
    }
    if(new_c->n_point>MAX_NUM_TRACING)
    {
      printf("The points of the curve has %zu\n", new_c->n_point);
      fprintf(stderr,"WARNING: please check the points number in the curve.\n");
      exit(EXIT_FAILURE);
    }
  }
  // write_curve("DEBUG_tracing_c", new_c);

  return new_c;
}

// tracing from a grid zone start point and finally arrive to the end curve.
// Here we assume the tracing direction HAS BEEN CORRECTED.
// We also asumme that the tracing line will arrive at the end curve.
// after use, do not forget to FREE it!!!!!
// We also asumme there is only ONE intersection on the end curve.
// THIS IS EXPECTED FOR 3D MAGNETIC FIELD. 
static Curve* create_2d_GridTubeCurve_3d_tracing(double start_R, double start_Z, double start_phi,Curve* end_curve, 
                                              ode_function* func,ode_solver* solver)
{
  if(func->ndim!=3)
  {
    fprintf(stderr, "create_2d_GridTubeCurve_2d_tracing used for 2D magnetic field.\n");
    exit(EXIT_FAILURE);
  }
  if (end_curve->n_point < 2) 
  {
    fprintf(stderr, "Error: end_curve has too few points.\n");
    exit(EXIT_FAILURE);
  }
  Curve* new_c=create_curve(MAX_NUM_TRACING);
  add_last_point_curve(new_c, start_R, start_Z);
  double curr_p[3]={start_R, start_Z, start_phi};
  double next_p[3];
  double t=0.00;
  double step_size=0.1;
  while(true)
  {
    solver->next_step(step_size, &t, curr_p, next_p, solver->solver_data, func);
    bool intersection=false;
    bool isclosed=false;
    for(int i=1; i<end_curve->n_point;i++)
    {
      // return 0 means found intersection
      // new_c->n_point>2 otherwise for the core region the first element of curve is always
      // has intersection with end curve.
      if(has_intersection(curr_p[0],curr_p[1],next_p[0],next_p[1],
                          end_curve->points[i-1].x, end_curve->points[i-1].y,
                          end_curve->points[i].x, end_curve->points[i].y)==0
                          && new_c->n_point>2) 
      {
        double intsect_x, intsect_y;
        get_intersection_point(curr_p[0],curr_p[1],next_p[0],next_p[1],
                               end_curve->points[i-1].x, end_curve->points[i-1].y,
                               end_curve->points[i].x, end_curve->points[i].y,
                               &intsect_x, &intsect_y);
        if(fabs(intsect_x-start_R)<1.0E-8&&fabs(intsect_y-start_Z)<1.0E-8)
        {
          add_last_point_curve(new_c,start_R, start_Z);
        }
        else
        {
          add_last_point_curve(new_c, intsect_x, intsect_y);
        }
        intersection=true;
        break;
      }
    }
    if(intersection||isclosed)
    {
      break;
    }
    else
    {
      add_last_point_curve(new_c, next_p[0], next_p[1]);
      curr_p[0]=next_p[0];
      curr_p[1]=next_p[1];
    }
    if(new_c->n_point>MAX_NUM_TRACING)
    {
      printf("The points of the curve has %zu\n", new_c->n_point);
      fprintf(stderr,"WARNING: please check the points number in the curve.\n");
      exit(EXIT_FAILURE);
    }
  }
  // write_curve("DEBUG_tracing_c", new_c);

  //Restore dim
  return new_c;
}

// the generation are grid is divided by generating the points along the magnetic surfaces one by one
// The corespoding structure is called grid tube. the fisrt line and grid point is known then calculate 
// the points along the second line of grid tube.
void generate_CARRE_2Dgrid_default(TwoDimGrid* grid,
                                   GridZone* gridzone,
                                   ode_function* func,
                                   ode_solver* solver)
{
  int np=gridzone->first_gridpoint_curve->n_point; //poloidal number along a magnetic field line.

  Curve* prev_c=gridzone->first_bnd_curve;
  Curve* prev_gpt_c=gridzone->first_gridpoint_curve;

  double* len_prev_gpt_c=malloc(np*sizeof(double));
  len_prev_gpt_c[0]=0.0;
  len_prev_gpt_c[np-1]=total_length_curve(prev_c);

  double d1=0.0;

  for(int i=1;i<np-1;i++)
  {
    d1 = ruban_curve(prev_c, &(prev_gpt_c->points[i]), d1);
    len_prev_gpt_c[i]=d1;
  }

  // write_array(len_prev_gpt_c, np, "DEBUG_len_prev_gpt_c");
  // write_curve("DEBUG_prev_gpt_c",prev_gpt_c);
  // write_curve("DEBUG_prev_c",prev_c);
  // ==== DO NOT FORGET CHECK and CORRECT the poloidal tracing direction
  check_poloidal_direction(gridzone, func, solver);

  Curve* curr_c=create_2d_GridTubeCurve_2d_tracing(gridzone->start_point_R[1],
                                                gridzone->start_point_Z[1],
                                                gridzone->end_curve,
                                                func, solver);
  Curve* curr_gpt_c=create_curve(np);

  expand_curve_size_with_NaN(curr_gpt_c,np);

  double* len_curr_gpt_c=malloc(np*sizeof(double));
  GirdTubeStr* gridtube=create_GridTube(prev_c, prev_gpt_c, len_prev_gpt_c, 
                                        curr_c, curr_gpt_c, len_curr_gpt_c,
                                        gridzone->guard_start[0],
                                        gridzone->guard_end[0],
                                        gridzone->pasmin[0]);

  // printf("%s %s %s\n","guard_start", "guard_end","pasmin");
  // printf("%lf %lf %lf\n", gridzone->guard_start[0], gridzone->guard_end[0],gridzone->pasmin[0]);

  /*****************************************************
  * CORE ALGORITHM: Calculate the poloidal distribution.
  *****************************************************/
  calc_points_from_CARRE(gridtube);

  char name[32];
  sprintf(name,"%s1",gridzone->name);
  write_curve(name,curr_gpt_c);

  for(int i=2;i<gridzone->nr; i++)
  {
  //Becareful!!!, Just change the adress and not copy the content.
    prev_c=curr_c;
    prev_gpt_c = curr_gpt_c;

    memcpy(len_prev_gpt_c,len_curr_gpt_c,np*sizeof(double));

    curr_c=create_2d_GridTubeCurve_2d_tracing(gridzone->start_point_R[i],
                                           gridzone->start_point_Z[i],
                                           gridzone->end_curve,
                                           func, solver);
    curr_gpt_c=create_curve(np);
    for(int j=0;j<np;j++)
    {
      add_last_point_curve(curr_gpt_c, 0.0, 0.0);
    }

    gridtube->prev_c=prev_c;
    gridtube->prev_gpt_c=prev_gpt_c;
    gridtube->len_prev_gpt_c=len_prev_gpt_c;
    gridtube->curr_c=curr_c;
    gridtube->curr_gpt_c=curr_gpt_c;
    gridtube->len_curr_gpt_c=len_curr_gpt_c;
    gridtube->guard_top=gridzone->guard_start[i];
    gridtube->guard_end=gridzone->guard_end[i];
    gridtube->pasmin = gridzone->pasmin[i];

  /*****************************************************
  * CORE ALGORITHM: Calculate the poloidal distribution.
  *****************************************************/
    calc_points_from_CARRE(gridtube);
  
  #ifdef DEBUG
    sprintf(name,"%s%d",gridzone->name, i);
    write_curve(name,curr_gpt_c);
  #endif
    
    free_curve(prev_c);
    free_curve(prev_gpt_c);
    if(i==gridzone->nr-1&&gridzone->sec_bnd)
    {
      //TODO specific operation for multiple X-points situations.
    }
  }

  //NOT FREE the entries in gridtube, but the gridtube space itself!
  free(gridtube);
  free_curve(curr_c);
  free_curve(curr_gpt_c);
  free(len_prev_gpt_c);
  free(len_curr_gpt_c);
  return;
}

//phim, nphi and array phi will decided the nfirst and nlast. 
//nstart indicate the firt nstart+1 number, nend indictae the last nend+1 number, are fixed positions.
static void calc_nfirst_nlast(double phim, int nphi, double* phi, int* nfirst_ptr, int* nlast_ptr)
{
  int idx_phim=-1;
  for(int i=0;i<nphi;i++)
  {
    if(fabs(phi[i]-phim)<EPS_TDGG)
    {
      idx_phim=i;
      break;
    }
  }
  if(idx_phim==-1)
  {
    fprintf(stderr,"phim %lf is not in the phi range\n",phim);
    exit(EXIT_FAILURE);
  }
  *nfirst_ptr = nphi - 1 - idx_phim;
  *nlast_ptr = idx_phim;
}

static void restore_3D_mag_direction(ode_function* func)
{
  if(func->ndim!=3)
  {
    fprintf(stderr,"restore_3D_mag_direction for 3D magnetic.\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0;i<3;i++)
  {
    func->rescale[i]=1.0;
  }
}

static void reverse_3D_mag_direction(ode_function* func)
{
  if(func->ndim!=3)
  {
    fprintf(stderr,"reverse_3D_mag_direction for 3D magnetic.\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0;i<3;i++)
  {
    func->rescale[i]=-1.0;
  }
}
// Fast but dangerous 3D line tracing.
// Start from pt, then transport along magnetic field line,
// until the pt[2] equals phi_tgt. Then store the value in pt_tgt.
// This function is very dangerous, the user should ensure the magnetic field is ok.
// 0 means arrive to the phi_tgt, 1 means failure.
// len_rz is used to record the total lengh in the RZ plane.
static int fast_3D_line_tracing(double* pt, double* pt_tgt, double phi_tgt, double* len_RZ,
                              ode_function* func, ode_solver* solver)
{
  *len_RZ=0.0;
  if(func->ndim!=3)
  {
    fprintf(stderr,"fast_3D_line_tracing only used for 3D magnetic field.\n");
    exit(EXIT_FAILURE);
  }
  int n=0;
  double pt_tmp[3];
  // Copy initial point
  for(int i=0;i<3;i++)
  {
    pt_tmp[i]=pt[i];
  }
  double next_pt_tmp[3];
  double t_tmp=0.0;
  
  while(true)
  {
    // Take one integration step
    solver->next_step(solver->step_size, &t_tmp, pt_tmp, next_pt_tmp, 
                      solver->solver_data, func);
    t_tmp=t_tmp + solver->step_size;
    n=n+1;
    *len_RZ=hypot(next_pt_tmp[0]-pt_tmp[0],next_pt_tmp[1]-pt_tmp[1])+*len_RZ;
    // Check if we reached the target phi value
    if(fabs(next_pt_tmp[2]-phi_tgt)<EPS_TDGG)
    {
      // #ifdef DEBUG
      //   printf("======================================================\n");
      //   printf("DEBUG Target          Phi: %.12f\n",phi_tgt);
      //   printf("DEBUG Start point R Z Phi: %.12f %.12f %.12f\n",pt[0], pt[1], pt[2]);
      //   printf("DEBUG End   point R Z Phi: %.12f %.12f %.12f\n",next_pt_tmp[0], next_pt_tmp[1], next_pt_tmp[2]);
      //   printf("======================================================\n");
      // #endif
      // Copy result to output array
      for(int i=0;i<3;i++)
      {
        pt_tgt[i]=next_pt_tmp[i];
      }
      return 0; // Success
    }

    // Update starting point for next iteration
    for(int j=0;j<3;j++)
    {
      pt_tmp[j]=next_pt_tmp[j];
    }

    // printf("DEBUG %12.f\n",pt_tmp[2]);

    // Check if we exceeded maximum number of steps
    if(n>MAX_NUM_TRACING)
    {
      fprintf(stderr,"Reached the tracing limitation %i\n",n);
      fprintf(stderr,"Please check all inputs for fast_3D_line_tracing\n");
      return 1; // Failure
    }
  }
}

// Fast but dangerous 2D line tracing.
// Start from pt, then transport along magnetic field line,
// until intersect with the Curve end_curve. 
// The intersection point is stored in intersct_pt
// This function is very dangerous, the user should ensure the magnetic field is ok.
// 0 means found the intersection, 1 means failure.
static int fast_3D_line_tracing_intersection(double* pt, Curve* end_curve, double* intersct_pt,
                                          ode_function* func, ode_solver* solver)
{
  int n=0;
  double curr_p[3];
  // Copy initial point
  for(int i=0;i<3;i++)
  {
    curr_p[i]=pt[i];
  }
  double next_p[3];
  double t_tmp=0.0;

  while (true)
  {
    solver->next_step(solver->step_size, &t_tmp, curr_p, next_p, solver->solver_data, func);
    n++;
    bool intersection=false;
    int np=end_curve->n_point;
    for(int i=1; i<np;i++)
    {
      // return 0 means found intersection
      if(has_intersection(curr_p[0],curr_p[1],next_p[0],next_p[1],
                          end_curve->points[i-1].x, end_curve->points[i-1].y,
                          end_curve->points[i].x, end_curve->points[i].y)==0) 
      {
        get_intersection_point(curr_p[0],curr_p[1],next_p[0],next_p[1],
                               end_curve->points[i-1].x, end_curve->points[i-1].y,
                               end_curve->points[i].x, end_curve->points[i].y,
                               &intersct_pt[0], &intersct_pt[1]);
        // Calculate interpolation ratio - position of intersection on curr_p to next_p line segment
        double dx = next_p[0] - curr_p[0];
        double dy = next_p[1] - curr_p[1];
        double line_length = sqrt(dx*dx + dy*dy);
    
        if(line_length > 1e-12) {  // Avoid division by zero
          double dx_to_intersect = intersct_pt[0] - curr_p[0];
          double dy_to_intersect = intersct_pt[1] - curr_p[1];
          double intersect_length = sqrt(dx_to_intersect*dx_to_intersect + dy_to_intersect*dy_to_intersect);
          double t_ratio = intersect_length / line_length;
          // Linear interpolation to calculate phi coordinate
          intersct_pt[2] = curr_p[2] + t_ratio * (next_p[2] - curr_p[2]);
        } else {
          intersct_pt[2] = curr_p[2];  // If step size is very small, use current point's Z coordinate
        }
        intersection=true;
        break;
      }
    }
    if(intersection)
    {
      return 0;
    }

    if(n>MAX_NUM_TRACING)
    {
      fprintf(stderr,"Unexpected: Reached the tracing limitation in fast_3D_line_tracing_intersection.\n");
      return 1;
    }
    for(int i=0;i<3;i++) 
    {
      curr_p[i] = next_p[i];
    }
  }
}


void update_sn_SepDistStr_PolSegmsInfo_EMC3_2Dgrid(PolSegmsInfo *polseg, SepDistStr* sepdist,
                                                   ode_function* func,ode_solver* solver,
                                                   double phim, int nphi, double* phi)

{
  if(!polseg||!sepdist)
  {
    fprintf(stderr,"Empty input for update_SepDist_PolSegmsInfo for EMC3 2Dgrid.\n");
    exit(EXIT_FAILURE);
  }
  if(strcmp(polseg->topo,"SNL")!=0)
  {
    fprintf(stderr,"Only support SNL topology.\n");
    exit(EXIT_FAILURE);
  }
  int nfirst;
  int nlast;
  calc_nfirst_nlast(phim, nphi, phi, &nfirst, &nlast);
  //!idx_phim is the index, nlast is the i-th
  int idx_phim=nlast;


  if(polseg->polsegments[0]->n_points<=nlast  //outer leg
     ||polseg->polsegments[1]->n_points<=nfirst) //inner leg
  {
    fprintf(stderr,"The resolution in phi direction is too high!\n");
    fprintf(stderr,"Decrease phi resolution (nphi and phi range) or increase PolSegStr resolution!\n");
    exit(EXIT_FAILURE);
  }

  if(polseg->polsegments[0]->n_points-nlast<10  //outer leg
     ||polseg->polsegments[1]->n_points-nfirst<10) //inner leg
  {
    printf("Warning: The resolution in phi direction is close to PolSegStr resolution.\n");
    printf("Suggest decrasing phi resolution and(or) increase PolSegStr resolution!\n");
  }

/**************************************************************************
* Check the magnetic field to ensure it can be used for 3D grid generation
***************************************************************************/
  //check the dimension
  if(func->ndim!=3)
  {
    fprintf(stderr,"The magnetic field is not 3D!\n");
    fprintf(stderr,"Even though is 2D grid generation, but it is the base of 3D grid.\n");
    fprintf(stderr,"Please check the magnetic field used for 3D grid genetration!\n");
    exit(EXIT_FAILURE);
  }

  //check the directions Br,Bz,Bphi whether they are consistent with direction definition.
  //restore rescale
  for(int i=0;i<3;i++)
  {
    func->rescale[i]=1.0;
  }
  //check Bphi(Toroidal)
  MagFieldTorSys *mag_field = (MagFieldTorSys *) func->data;
  if(mag_field->b0r0>0)
  {
    fprintf(stderr,"BT is into the page (away from you).\n");
    fprintf(stderr,"The direction of BT is not support.\n");
    exit(EXIT_FAILURE);
  }
  //start point for checking the Br and Bz.
  //We only cheeck at the inner target and assume that if inner is ok than outer is ok.
  double p1[3];
  int idx=sepdist->index[0];
  int n_point=sepdist->edges[idx]->gridpoint_curve->n_point;
  p1[0]=sepdist->edges[idx]->gridpoint_curve->points[n_point-1].x;
  p1[1]=sepdist->edges[idx]->gridpoint_curve->points[n_point-1].y;
  p1[2]=phi[nphi-1];

  double t_tmp=0.0;
  double p2[3];
  solver->next_step(solver->step_size, &t_tmp, p1, p2, solver->solver_data, func);
  printf("Start Point                  R Z Phi %.12f %.12f %.12f\n",p1[0],p1[1],p1[2]);
  printf("After One step tracing Point R Z Phi %.12f %.12f %.12f\n",p2[0],p2[1],p2[2]);
  double p3[2];
  double p4[2];
  p3[0]=sepdist->edges[idx]->gridpoint_curve->points[n_point-1].x;
  p3[1]=sepdist->edges[idx]->gridpoint_curve->points[n_point-1].y;
  p4[0]=sepdist->edges[idx]->gridpoint_curve->points[n_point-1-1].x;
  p4[1]=sepdist->edges[idx]->gridpoint_curve->points[n_point-1-1].y;
  if(dot_product(p2[0]-p1[0], p2[1]-p1[1], p4[0]-p3[0], p4[1]-p3[1])<0.0)
  {
    fprintf(stderr,"Ip is out of the page (towards you).\n");
    fprintf(stderr,"Thus, Br&Bz is not consistent with expected direction.\n");
    fprintf(stderr,"The direction of Ip is not support.\n");
    exit(EXIT_FAILURE);
  }
  printf("BT and Ip directions are consitent with direction definition.\n");

/***********************************************
*  Correct the inner leg ((sepdist->edges[index[0]])  *
************************************************/

  int idx_tmp=sepdist->index[0];

  if(sepdist->edges[idx_tmp]->gridpoint_curve->n_point!=sepdist->edges[idx_tmp]->n_size)
  {
    fprintf(stderr,"The number of points of gridcurve is not consistent with size of normal distribution.\n");
    exit(EXIT_FAILURE);
  }
  
  //new gridpoint_curve for inner leg will used to replace the current one;
  //new normalized distribution for inner leg will be updated;

  int npoint_tmp=sepdist->edges[idx_tmp]->gridpoint_curve->n_point;
  // DO NOT FORGET FREE/REPLACE
  Curve* new_gpc_inner=create_curve(npoint_tmp);
  expand_curve_size_with_NaN(new_gpc_inner, npoint_tmp);

  // Calculate the new gridpoint_curve [npoint_tmp-1]
  set_point_curve(new_gpc_inner,npoint_tmp-1, 
                  sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].x,
                  sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].y);
  
  //BECAREFUL, nfirst is veiw from inner target to outer target along magnetic filed.
  //The direction of gridpointcurve is from X-point to inner/outer target.
  //So, the nfirst+1 points for inner leg of the last points of gridpointcurve.
  //the index of the last points are [npoint_tmp-nfirst-1:npoint_tmp-1].
  //Calculate the fixed point from the last [npoint_tmp-nfirst-1] to the last [npoint_tmp-2].
  //The last point index is [npoint_tmp-1] and no need to change.


  // Calculate the new gridpoint_curve from [npoint_tmp-nfirst-1:npoint_tmp-2]
  restore_3D_mag_direction(func);
  for(int i=1; i<nfirst+1; i++)
  {
    double pt[3];
    double pt_tgt[3];
    pt[0]=sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].x;
    pt[1]=sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].y;
    pt[2]=phi[idx_phim+i];
    printf("DEBUG Phi is %.12f\n",pt[2]);
    double len_RZ;
    fast_3D_line_tracing(pt,pt_tgt, phim,&len_RZ, func, solver);
    set_point_curve(new_gpc_inner,npoint_tmp-1-i, pt_tgt[0],pt_tgt[1]);
  }

  printf("DEBUG %.12f %.12f\n",new_gpc_inner->points[npoint_tmp-nfirst-1].x,new_gpc_inner->points[npoint_tmp-nfirst-1].y);

  // Create temporary new separatrix line DLList.
  DLListNode* head_in_tmp=copy_DLList(sepdist->edges[idx_tmp]->head);
  if(insert_point_for_DLList(head_in_tmp, 
                             new_gpc_inner->points[npoint_tmp-nfirst-1].x,
                             new_gpc_inner->points[npoint_tmp-nfirst-1].y))
  {
    fprintf(stderr,"The point in not on the separatrix line.\n");
    exit(EXIT_FAILURE);
  }
  write_DLList(head_in_tmp,"DEBUG_SEPLINE_IN");

  cut_DLList_after_point(head_in_tmp,
                                new_gpc_inner->points[npoint_tmp-nfirst-1].x,
                                new_gpc_inner->points[npoint_tmp-nfirst-1].y);

  // Create the coresponding normal distribution for the new separatrix line DLList
  double* normdist_in_tmp=malloc((npoint_tmp-nfirst)*sizeof(double));
  double norm_factor = sepdist->edges[idx_tmp]->norm_dist[npoint_tmp-nfirst-1];

  if(fabs(norm_factor) < EPS_TDGG) 
  {
    fprintf(stderr,"Normalization factor is too small or zero!\n");
    exit(EXIT_FAILURE);
  }  
  for(int i=0;i<npoint_tmp-nfirst;i++)
  {
    normdist_in_tmp[i]=sepdist->edges[idx_tmp]->norm_dist[i]/norm_factor;
    printf("DEBUG normdist_in_tmp %i %lf\n",i,normdist_in_tmp[i]);
  }
  
  // Calculate the new gridpoint_curve from [0:npoint_tmp-nfirst-1]
  Curve* gpc_in_tmp=create_gridpoint_curve(head_in_tmp, normdist_in_tmp, npoint_tmp-nfirst);
  write_curve("DEBUG_gpc_in_tmp",gpc_in_tmp);

  //fill out the new_gpc_inner from gpc_in_tmp[0:npoint_tmp-nfirst-2]
  for(int i=0; i<npoint_tmp-nfirst-1; i++)
  {
    new_gpc_inner->points[i].x=gpc_in_tmp->points[i].x;
    new_gpc_inner->points[i].y=gpc_in_tmp->points[i].y;
  }

  //Update the normalized distribution,norm_dist[0]=0.0 and norm_dist[npoint-1]=1.0
  double tot_len_tmp = total_length_curve(new_gpc_inner);
  for(int i=0; i<npoint_tmp;i++)
  {
    sepdist->edges[idx_tmp]->norm_dist[i]=length_curve(new_gpc_inner,i+1)/tot_len_tmp;
    //the 2nd polseg is coresponding to the inner leg!!! 
    polseg->polsegments[1]->norm_dist[i]=sepdist->edges[idx_tmp]->norm_dist[i];
    printf("DEBUG new norm dist %i %.12f\n",i, sepdist->edges[idx_tmp]->norm_dist[i]);
  }
  polseg->polsegments[1]->norm_dist[0]=0.0;
  polseg->polsegments[1]->norm_dist[npoint_tmp-1]=1.0;

  //Replace the old curve by the new curve.
  free_curve(sepdist->edges[idx_tmp]->gridpoint_curve);
  sepdist->edges[idx_tmp]->gridpoint_curve=new_gpc_inner;

  write_curve("DEBUG_new_gpc_in",new_gpc_inner);
  free_curve(gpc_in_tmp);
  free_DLList(head_in_tmp);
  free(normdist_in_tmp);

/*****************************************************
*  Correct the outer leg (sepdist->edges[index[1]])   *
*****************************************************/
  idx_tmp=sepdist->index[1];

  if(sepdist->edges[idx_tmp]->gridpoint_curve->n_point!=sepdist->edges[idx_tmp]->n_size)
  {
    fprintf(stderr,"The number of points of gridcurve is not consistent with size of normal distribution.\n");
    exit(EXIT_FAILURE);
  }

  //new gridpoint_curve for outer leg will used to replace the current one;
  //new normalized distribution for outer leg will be updated;
  npoint_tmp=sepdist->edges[idx_tmp]->gridpoint_curve->n_point;
  // DO NOT FORGET FREE/REPLACE
  Curve* new_gpc_outer=create_curve(npoint_tmp);
  expand_curve_size_with_NaN(new_gpc_outer, npoint_tmp);

  // Calculate the new gridpoint_curve [npoint_tmp-1]
  set_point_curve(new_gpc_outer,npoint_tmp-1, 
                  sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].x,
                  sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].y);

  //BECAREFUL, nlast is veiw from inner target to outer target along magnetic filed.
  //The direction of gridpointcurve is from X-point to inner/outer target.
  //So, the nlast+1 points for outer leg of the last points of gridpointcurve.
  //the index of the last points are [npoint_tmp-nlast-1:npoint_tmp-1].
  //Calculate the fixed point from the last [npoint_tmp-nlast-1] to the last [npoint_tmp-2].
  //The last point index is [npoint_tmp-1] and no need to change.


  //For outer leg, the tracing direction is along the reversed magnetic field line.
  reverse_3D_mag_direction(func);
  // Calculate the new gridpoint_curve from [npoint_tmp-nlast1-1:npoint_tmp-2]
  for(int i=1; i<nlast+1; i++)
  {
    double pt[3];
    double pt_tgt[3];
    pt[0]=sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].x;
    pt[1]=sepdist->edges[idx_tmp]->gridpoint_curve->points[npoint_tmp-1].y;
    pt[2]=phi[idx_phim-i];
    printf("DEBUG Phi is %.12f\n",pt[2]);
    double len_RZ;
    fast_3D_line_tracing(pt,pt_tgt, phim,&len_RZ, func, solver);
    set_point_curve(new_gpc_outer,npoint_tmp-1-i, pt_tgt[0],pt_tgt[1]);
  }
  restore_3D_mag_direction(func);

  printf("DEBUG %.12f %.12f\n",
          new_gpc_outer->points[npoint_tmp-nlast-1].x,
          new_gpc_outer->points[npoint_tmp-nlast-1].y);

  // Create temporary new separatrix line DLList.
  DLListNode* head_out_tmp=copy_DLList(sepdist->edges[idx_tmp]->head);
  if(insert_point_for_DLList(head_out_tmp, 
                             new_gpc_outer->points[npoint_tmp-nlast-1].x,
                             new_gpc_outer->points[npoint_tmp-nlast-1].y))
  {
    fprintf(stderr,"The point in not on the separatrix line.\n");
    exit(EXIT_FAILURE);
  }

  cut_DLList_after_point(head_out_tmp,
                               new_gpc_outer->points[npoint_tmp-nlast-1].x,
                               new_gpc_outer->points[npoint_tmp-nlast-1].y);
  write_DLList(head_out_tmp,"DEBUG_SEPLINE_OUT");


  // Create the coresponding normal distribution for the new separatrix line DLList
  double* normdist_out_tmp=malloc((npoint_tmp-nlast)*sizeof(double));
  norm_factor = sepdist->edges[idx_tmp]->norm_dist[npoint_tmp-nlast-1];

  if(fabs(norm_factor) < EPS_TDGG) 
  {
    fprintf(stderr,"Normalization factor is too small or zero!\n");
    exit(EXIT_FAILURE);
  }  

  for(int i=0;i<npoint_tmp-nlast;i++)
  {
    normdist_out_tmp[i]=sepdist->edges[idx_tmp]->norm_dist[i]/norm_factor;
    printf("DEBUG normdist_out_tmp %i %.12f\n",i,normdist_out_tmp[i]);
  }

  // Calculate the new gridpoint_curve from [0:npoint_tmp-nfirst-1]
  Curve* gpc_out_tmp=create_gridpoint_curve(head_out_tmp, normdist_out_tmp, npoint_tmp-nlast);
  write_curve("DEBUG_gpc_out_tmp",gpc_out_tmp);

  //fill out the new_gpc_outer from gpc_out_tmp[0:npoint_tmp-nfirst-2]
  for(int i=0; i<npoint_tmp-nlast-1; i++)
  {
    new_gpc_outer->points[i].x=gpc_out_tmp->points[i].x;
    new_gpc_outer->points[i].y=gpc_out_tmp->points[i].y;
  }

  //Update the normalized distribution,norm_dist[0]=0.0 and norm_dist[npoint-1]=1.0
  tot_len_tmp = total_length_curve(new_gpc_outer);
  for(int i=0; i<npoint_tmp;i++)
  {
    sepdist->edges[idx_tmp]->norm_dist[i]=length_curve(new_gpc_outer,i+1)/tot_len_tmp;
    //the 1st polseg is coresponding to the inner leg!!! 
    polseg->polsegments[0]->norm_dist[i]=sepdist->edges[idx_tmp]->norm_dist[i];
    printf("DEBUG new norm dist %i %.12f\n",i, sepdist->edges[idx_tmp]->norm_dist[i]);
  }
  polseg->polsegments[0]->norm_dist[0]=0.0;
  polseg->polsegments[0]->norm_dist[npoint_tmp-1]=1.0;

  //Replace the old curve by the new curve.
  free_curve(sepdist->edges[idx_tmp]->gridpoint_curve);
  sepdist->edges[idx_tmp]->gridpoint_curve=new_gpc_outer;

  write_curve("DEBUG_new_gpc_out",new_gpc_outer);
  free_curve(gpc_out_tmp);
  free_DLList(head_out_tmp);
  free(normdist_out_tmp);

}


void generate_EMC3_2Dgrid_default(TwoDimGrid* grid,
                                   GridZone* gridzone,
                                   ode_function* func,
                                   ode_solver* solver,
                                   double phim, int nphi, double* phi)
{
/*****************************************************
*  Check input variables
*****************************************************/
  if(!phi||nphi<2)
  {
    fprintf(stderr,"The inputs about phi are wrong.\n");
    exit(EXIT_FAILURE);
  }
  int nfirst, nlast, idx_phim;
  // For core gridzone
  if(strncmp(gridzone->name,"CORE",4)==0)
  {
    nfirst=0;
    nlast=0;
    idx_phim=0;
  }
  // For SOL and PFR gridzone
  else
  {
    calc_nfirst_nlast(phim, nphi, phi, &nfirst, &nlast);
    idx_phim=nlast;
    #ifdef DEBUG
      printf("nfirst %i, nlast %i\n",nfirst,nlast);
    #endif
  }
  
/**************************************************************************
* Check the magnetic field to ensure it can be used for 3D grid generation
***************************************************************************/
 //TO DO

/**************************************************************************
* New first boundary curve and new first gridpoint curve
* This is because the first nfirst+1 and the (nlast+1) to the last
***************************************************************************/

  int np=gridzone->first_gridpoint_curve->n_point-nfirst-nlast; //poloidal number along a magnetic field line.
  
  //create a temperary dllist for first_bnd_curve 
  //then instert fixed grid point point[nfirst] and point[nlast].
  //Becareful grid point point[nfirst] is the (nfirst+1)-th, point[nlast10] is the (nlast+1)-th TO THE LAST
  DLListNode* fbc_tmp_dllist = create_DLListNode(gridzone->first_bnd_curve->points[0].x,
                                                 gridzone->first_bnd_curve->points[0].y);
  DLListNode* tail_tmp = fbc_tmp_dllist;
  int npoint_fbc=gridzone->first_bnd_curve->n_point;
  for(int i=1;i<npoint_fbc;i++)
  {
    add_DLListnode_at_tail(&tail_tmp, 
                            gridzone->first_bnd_curve->points[i].x,
                            gridzone->first_bnd_curve->points[i].y);
  }

  if(insert_point_for_DLList(fbc_tmp_dllist, gridzone->first_gridpoint_curve->points[nfirst].x, 
                                             gridzone->first_gridpoint_curve->points[nfirst].y))
  {
    fprintf(stderr,"The grid point (nfirst+1)-th  in not on the first_bnd_curve.\n");
    exit(EXIT_FAILURE);
  };
  cut_DLList_before_point(&fbc_tmp_dllist, gridzone->first_gridpoint_curve->points[nfirst].x, 
                                           gridzone->first_gridpoint_curve->points[nfirst].y);


  int tmp = gridzone->first_gridpoint_curve->n_point-1-nlast;
  if(insert_point_for_DLList(fbc_tmp_dllist, gridzone->first_gridpoint_curve->points[tmp].x, 
                                             gridzone->first_gridpoint_curve->points[tmp].y))
  {
    fprintf(stderr,"The grid point (nfirst+1)-th to the last in not on the first_bnd_curve.\n");
    exit(EXIT_FAILURE);
  };
  cut_DLList_after_point(fbc_tmp_dllist, gridzone->first_gridpoint_curve->points[tmp].x, 
                                         gridzone->first_gridpoint_curve->points[tmp].y);

  #ifdef DEBUG
    write_DLList(fbc_tmp_dllist,"fbc_tmp_dllist");
  #endif
  Curve* new_first_bnd_curve=create_curve(npoint_fbc);
  DLListNode* head_tmp=fbc_tmp_dllist;
  while(head_tmp)
  {
    add_last_point_curve(new_first_bnd_curve, head_tmp->r,head_tmp->z);
    head_tmp = head_tmp->next;
  }
  Curve* new_first_gridpoint_curve=create_curve(np);

  tmp=gridzone->first_gridpoint_curve->n_point-nlast;
  for(int i=nfirst;i<tmp;i++)
  {
    add_last_point_curve(new_first_gridpoint_curve,
                         gridzone->first_gridpoint_curve->points[i].x,
                         gridzone->first_gridpoint_curve->points[i].y);
  }

  #ifdef DEBUG
    write_curve("new_first_bnd_curve",new_first_bnd_curve);
    write_curve("new_first_gridpoint_curve",new_first_gridpoint_curve);
  #endif

  free_DLList(fbc_tmp_dllist);

/**************************************************************************
* New senond boundary curve and new second gridpoint curve
* This is because the first nfirst+1 and the last 
***************************************************************************/
//TODO
  if(gridzone->sec_bnd)
  {
    //Todo
    fprintf(stderr,"Not yet second boundary curve.\n");
    exit(EXIT_FAILURE);
  }

/**************************************************************************
* New start_point_R, start_point_Z, guard_start,
* This is because the first (nfirst+1) points 
***************************************************************************/
  int nr =gridzone->nr;
  double* new_start_point_R=malloc(nr*sizeof(double));
  double* new_start_point_Z=malloc(nr*sizeof(double));
  double* new_guard_start=malloc(nr*sizeof(double));
  double* new_guard_end=malloc(nr*sizeof(double));

  //to store the grid points decided by line tracing. The so-called first part gridpoint
  double*** first_gridpoint=allocate_3d_array(nfirst+1,nr,2);

  //restore the diection
  for(int i=0;i<3;i++)
  {
    func->rescale[i]=1.0;
  }
  
  //calculate first_gridpoint 
  for(int i=0;i<nr;i++)
  {
    double pt_tmp[3];
    double next_pt_tmp[3];
    double R_tmp = gridzone->start_point_R[i];
    double Z_tmp = gridzone->start_point_Z[i];
    double phi_tgt = phi[idx_phim];
    double len_RZ_tmp=0.0;
    for(int j=0;j<nfirst+1;j++)
    {
      if(j==0)
      {
        first_gridpoint[j][i][0]=R_tmp;
        first_gridpoint[j][i][1]=Z_tmp;
      }
      else
      {
        pt_tmp[0]=R_tmp;
        pt_tmp[1]=Z_tmp;
        pt_tmp[2]=phi[idx_phim+j];

        if(fast_3D_line_tracing(pt_tmp, next_pt_tmp, phi_tgt, &len_RZ_tmp, func, solver))
        {
          fprintf(stderr,"UNEXPECTED ERROR: Cannot Found the point in nr %i np %i\n", i, j);
          exit(EXIT_FAILURE);
        }
        first_gridpoint[j][i][0]=next_pt_tmp[0];
        first_gridpoint[j][i][1]=next_pt_tmp[1];
        //Update the new guard_start
        if(j==nfirst)
        {
          if(gridzone->guard_start[i]>len_RZ_tmp)
          {
            new_guard_start[i]=gridzone->guard_start[i]-len_RZ_tmp;
          }
          else
          {
            new_guard_start[i]=0.0;
          }
          #ifdef DEBUG
            printf("guard_start %.12f len_RZ %.12f new_guard_start %.12f\n",
                    gridzone->guard_start[i],len_RZ_tmp,new_guard_start[i]);
          #endif
        }
      }
    }
  }
  //consistency check.
  for(int i=0;i<nfirst+1;i++)
  {
    if(fabs(first_gridpoint[i][0][0]- gridzone->first_gridpoint_curve->points[i].x)>EPS_TDGG||
       fabs(first_gridpoint[i][0][1]- gridzone->first_gridpoint_curve->points[i].y)>EPS_TDGG)
    {
      fprintf(stderr,"UNEXPECTED ERROR: first gridpoint not consistent with gridpoint_curve\n");
      exit(EXIT_FAILURE);
    }
  }

  //update start_point_R and start_point_Z
  for(int i=0;i<nr;i++)
  {
    new_start_point_R[i]=first_gridpoint[nfirst][i][0];
    new_start_point_Z[i]=first_gridpoint[nfirst][i][1];
  }

  char name_tmp[32];
  sprintf(name_tmp,"%s_FIRST",gridzone->name);
  write_3d_array(nfirst+1,nr,2,first_gridpoint,name_tmp,2);

/**************************************************************************
* New guard_end and endcurve
* This is because the  (nlast1+1) points to the last 
***************************************************************************/
  //new end_curve is NOT the old end cuvre by tracing.
  //new end_curve is by the magnetic line intersection and 
  // then tracing from intersection to proper position.

  int n_new_endcurve=gridzone->nr;
  Curve* new_end_curve=create_curve(n_new_endcurve);

  
  //to store the grid points decided by line tracing. The so-called last part gridpoint
  double*** last_gridpoint=allocate_3d_array(nlast+1,nr,2);

  tmp=gridzone->first_bnd_curve->n_point;
  last_gridpoint[nlast][0][0]=gridzone->first_bnd_curve->points[tmp-1].x;
  last_gridpoint[nlast][0][1]=gridzone->first_bnd_curve->points[tmp-1].y;


  for(int i=1;i<nr;i++)
  {
    double pt_tmp[3];
    double phi_tgt = phi[idx_phim];
    pt_tmp[0]=new_start_point_R[i];
    pt_tmp[1]=new_start_point_Z[i];
    pt_tmp[2]=phi[idx_phim];
    double intersct_pt[3];
    if(fast_3D_line_tracing_intersection(pt_tmp, gridzone->end_curve, intersct_pt, func, solver))
    {
      fprintf(stderr,"UNEXPECTED ERROR: Cannot Found the intersecton.\n");
      exit(EXIT_FAILURE);
    }
    last_gridpoint[nlast][i][0]=intersct_pt[0];
    last_gridpoint[nlast][i][1]=intersct_pt[1];
  }

  //reversed magnetic field
  reverse_3D_mag_direction(func);
  for(int i=0;i<nr;i++)
  {
    double pt_tmp[3];
    double next_pt_tmp[3];
    double phi_tgt = phi[idx_phim];
    double len_RZ_tmp=0.0;

    for(int j=0;j<nlast;j++)
    {
      pt_tmp[0]=last_gridpoint[nlast][i][0];
      pt_tmp[1]=last_gridpoint[nlast][i][1];
      pt_tmp[2]=phi[j];
      // #ifdef DEBUG
      //   printf("pt_tmp R Z Phi: %.12f %.12f %.12f\n",pt_tmp[0],pt_tmp[1],pt_tmp[2]);
      //   printf("phi_tgt: %.12f\n",phi_tgt);
      // #endif
      if(fast_3D_line_tracing(pt_tmp, next_pt_tmp, phi_tgt, &len_RZ_tmp, func, solver))
      {
        fprintf(stderr,"UNEXPECTED ERROR: Cannot Found the point in nr %i np %i\n", i, j);
        exit(EXIT_FAILURE);
      }
      last_gridpoint[j][i][0]=next_pt_tmp[0];
      last_gridpoint[j][i][1]=next_pt_tmp[1];
      if(j==0)
      {
        if(gridzone->guard_end[i]>len_RZ_tmp)
        {
          new_guard_end[i]=gridzone->guard_end[i]-len_RZ_tmp;
        }
        else
        {
          new_guard_end[i]=0.0;
        }
        #ifdef DEBUG
            printf("guard_end %.12f len_RZ %.12f new_guard_end %.12f\n",
                    gridzone->guard_end[i],len_RZ_tmp,new_guard_end[i]);
        #endif
      }
    }
  }

  //restore the magnetic field
  restore_3D_mag_direction(func);

  // //consistency check. 
  tmp=gridzone->first_gridpoint_curve->n_point-1-nlast;
  sprintf(name_tmp,"%s_FIRST_GPC",gridzone->name);
  write_curve(name_tmp,gridzone->first_gridpoint_curve);
  for(int i=0;i<nlast+1;i++)
  {
    if(fabs(last_gridpoint[i][0][0]- gridzone->first_gridpoint_curve->points[tmp+i].x)>EPS_TDGG||
       fabs(last_gridpoint[i][0][1]- gridzone->first_gridpoint_curve->points[tmp+i].y)>EPS_TDGG)
    {
      fprintf(stderr,"UNEXPECTED ERROR: last gridpoint not consistent with gridpoint_curve\n");
      exit(EXIT_FAILURE);
    }
  }
  
  //Update new_end_curve
  for(int i=0;i<n_new_endcurve;i++)
  {
    add_last_point_curve(new_end_curve, 
                        last_gridpoint[0][i][0],
                        last_gridpoint[0][i][1]);
  }

  sprintf(name_tmp,"%s_new_end_curve",gridzone->name);
  write_curve(name_tmp, new_end_curve);

  sprintf(name_tmp,"%s_LAST",gridzone->name);
  write_3d_array(nlast+1,nr,2,last_gridpoint,name_tmp,2);


/**************************************************************************
* Firtst GridTube 
***************************************************************************/
// Similar to generate_CARRE_2Dgrid_default function

  Curve* prev_c=new_first_bnd_curve;
  Curve* prev_gpt_c=new_first_gridpoint_curve;
  double* len_prev_gpt_c=malloc(np*sizeof(double));
  len_prev_gpt_c[0]=0.0;
  len_prev_gpt_c[np-1]=total_length_curve(prev_c);

  double d1=0.0;

  for(int i=1;i<np-1;i++)
  {
    d1 = ruban_curve(prev_c, &(prev_gpt_c->points[i]), d1);
    len_prev_gpt_c[i]=d1;
  }
  double new_start_point_phi=phi[nlast];
  Curve* curr_c=create_2d_GridTubeCurve_3d_tracing(new_start_point_R[1],
                                                   new_start_point_Z[1],
                                                   new_start_point_phi,
                                                   new_end_curve,
                                                   func, solver);
  sprintf(name_tmp,"%s_curr_c",gridzone->name);
  write_curve(name_tmp,curr_c); 

  Curve* curr_gpt_c=create_curve(np);
  expand_curve_size_with_NaN(curr_gpt_c,np);
  
  double* len_curr_gpt_c=malloc(np*sizeof(double));
  GirdTubeStr* gridtube=create_GridTube(prev_c, prev_gpt_c, len_prev_gpt_c, 
                                        curr_c, curr_gpt_c, len_curr_gpt_c,
                                        new_guard_start[0],
                                        new_guard_end[0],
                                        gridzone->pasmin[0]);

  /*****************************************************
  * CORE ALGORITHM: Calculate the poloidal distribution.
  *****************************************************/
  calc_points_from_CARRE(gridtube);

  sprintf(name_tmp,"%s_2DBASE_1",gridzone->name);
  write_curve(name_tmp,curr_gpt_c);

  for(int i=2;i<gridzone->nr; i++)
  {
  //Becareful!!!, Just change the adress and not copy the content.
    prev_c=curr_c;
    prev_gpt_c = curr_gpt_c;

    memcpy(len_prev_gpt_c,len_curr_gpt_c,np*sizeof(double));

    curr_c=create_2d_GridTubeCurve_3d_tracing(new_start_point_R[i],
                                              new_start_point_Z[i],
                                              new_start_point_phi,
                                              new_end_curve,
                                              func, solver);
    curr_gpt_c=create_curve(np);
    for(int j=0;j<np;j++)
    {
      add_last_point_curve(curr_gpt_c, 0.0, 0.0);
    }

    gridtube->prev_c=prev_c;
    gridtube->prev_gpt_c=prev_gpt_c;
    gridtube->len_prev_gpt_c=len_prev_gpt_c;
    gridtube->curr_c=curr_c;
    gridtube->curr_gpt_c=curr_gpt_c;
    gridtube->len_curr_gpt_c=len_curr_gpt_c;
    gridtube->guard_top=new_guard_start[i];
    gridtube->guard_end=new_guard_end[i];
    gridtube->pasmin = gridzone->pasmin[i];

  /*****************************************************
  * CORE ALGORITHM: Calculate the poloidal distribution.
  *****************************************************/
    calc_points_from_CARRE(gridtube);
  
  #ifdef DEBUG
    sprintf(name_tmp,"%s_2DBASE_%d",gridzone->name, i);
    write_curve(name_tmp,curr_gpt_c);
  #endif
    
    free_curve(prev_c);
    free_curve(prev_gpt_c);

    if(i==gridzone->nr-1&&gridzone->sec_bnd)
    {
      //TODO specific operation for multiple X-points situations.
    }
  }

/**************************************************************************
* Free
***************************************************************************/
  free_curve(new_first_bnd_curve);
  free_curve(new_first_gridpoint_curve);
  free_curve(new_end_curve); 
  free(new_start_point_R);
  free(new_start_point_Z);
  free(new_guard_start);
  free(new_guard_end);
  free_3d_array(first_gridpoint);
  free_3d_array(last_gridpoint);

}