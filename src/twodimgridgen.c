#include "twodimgridgen.h"
#include "curve.h"
#include <math.h>
#include "mathbase.h"
#include "utils.h"

#define DEFAULT_MARGIN 20
#define MAX_NUM_TRACING 40000

#ifndef EPS_TDGG
#define EPS_TDGG 1.0E-10
#endif

#ifndef NRELAX
#define NRELAX 1000
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

TwoDimGrid* create_2Dgrid_default(int npol, int nrad)
{
    int margin_pol = DEFAULT_MARGIN;
    int margin_rad = DEFAULT_MARGIN;

    /* Allocate the grid container */
    TwoDimGrid* grid = malloc(sizeof(TwoDimGrid));
    if (!grid) 
    {
        fprintf(stderr, "Failed to allocate TwoDimGrid.\n");
        exit(EXIT_FAILURE);
    }

    /* Logical size */
    grid->npol = npol;
    grid->nrad = nrad;

    /* Physical capacity (size + padding on both sides) */
    grid->cap_npol  = npol + 2 * margin_pol;
    grid->cap_nrad  = nrad + 2 * margin_rad;

    /* Offsets that map (ir,ip) → physical index */
    grid->offset_pol = margin_pol;
    grid->offset_rad = margin_rad;

    /* Allocate and zero-initialize the point array */
    size_t total = (size_t)grid->cap_npol * grid->cap_nrad;
    grid->points = malloc(total * sizeof(GridPoint));
    if (!grid->points) 
    {
      fprintf(stderr, "Failed to allocate grid points.\n");
      free(grid);
      exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < total; ++i) 
    {
      grid->points[i].x = NAN;
      grid->points[i].y = NAN;
    }
    if (!grid->points) {
        fprintf(stderr, "Failed to allocate grid points.\n");
        free(grid);
        exit(EXIT_FAILURE);
    }
    return grid;
}

/*--------------------------------------------------------------------
  Release all memory held by a 2-D grid
--------------------------------------------------------------------*/
void free_2Dgrid(TwoDimGrid* grid)
{
    if (!grid) return;
    free(grid->points);
    free(grid);
}

/*--------------------------------------------------------------------
  Internal helper: convert logical (ir, ip) to a linear index
--------------------------------------------------------------------*/
static inline size_t idx_2Dgrid(const TwoDimGrid* g, int ir, int ip)
{
    return (size_t)(g->offset_rad + ir) * g->cap_npol
         + (size_t)(g->offset_pol + ip);
}

/*--------------------------------------------------------------------
  Return a pointer to the GridPoint at (ir, ip)
--------------------------------------------------------------------*/
GridPoint* get_point_2Dgrid(TwoDimGrid* g, int ir, int ip)
{
    return &g->points[idx_2Dgrid(g, ir, ip)];
}

/*--------------------------------------------------------------------
  Read-only accessors for x and y
--------------------------------------------------------------------*/
double get_x_2Dgrid(const TwoDimGrid* g, int ir, int ip)
{
    return g->points[idx_2Dgrid(g, ir, ip)].x;
}

double get_y_2Dgrid(const TwoDimGrid* g, int ir, int ip)
{
    return g->points[idx_2Dgrid(g, ir, ip)].y;
}

/*--------------------------------------------------------------------
  Set the coordinates of point (ir, ip)
--------------------------------------------------------------------*/
void set_point_2Dgrid(TwoDimGrid* g, int ir, int ip, double x, double y)
{
    GridPoint* p = &g->points[idx_2Dgrid(g, ir, ip)];
    p->x = x;
    p->y = y;
}



/****************************************************************************
*    INHERIT FROM CAREE
****************************************************************************/

static double calc_length(CurvePoint* p1, CurvePoint* p2) 
{
    return hypot(p1->x-p2->x, p1->y-p2->y);
}

static double dot_product(double x1, double y1, double x2, double y2) 
{
    return x1*x2+y1*y2;
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
    for(int i=0;i<n_point;i++)
    {
      add_last_point_curve(tmp_gpt_c,0.0,0.0);
    }

    tmp_length_points[0] = 0.0; 
    double length = total_length_curve(tube->curr_c);
    tmp_length_points[n_point-1] = length;

    double prev_length = tube->len_prev_gpt_c[n_point-1];
    
    // the first and the last point of mesh points are the same with curve
    int size = tube->curr_c->n_point-1;
    tube->curr_gpt_c->points[0].x=tube->curr_c->points[0].x;
    tube->curr_gpt_c->points[0].y=tube->curr_c->points[0].y;
    tube->curr_gpt_c->points[n_point-1].x=tube->curr_c->points[size-1].x;
    tube->curr_gpt_c->points[n_point-1].y=tube->curr_c->points[size-1].y;

    tmp_gpt_c->points[0].x = tube->curr_c->points[0].x;
    tmp_gpt_c->points[0].y = tube->curr_c->points[0].y;
    tmp_gpt_c->points[n_point-1].x = tube->curr_c->points[size-1].x;
    tmp_gpt_c->points[n_point-1].y = tube->curr_c->points[size-1].y;

    //DEBUG
    // write_curve("DEBUG_prev_gpt_c", tube->prev_gpt_c);
    // write_curve("DEBUG_prev_c",tube->prev_c);
    printf("DEBUG length: %lf\n",length);
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
          //printf("del: %.10f\n", del);
          //printf("ipol: %d, length_point %.10f\n",ipol, tube->length_points[ipol]);
          coordnates_in_curve(tube->curr_c, tube->len_curr_gpt_c[ipol], &(tube->curr_gpt_c->points[ipol]));

          ortmax=max(ortmax,fabs(ortho_current));
        }
        // printf("debug: ortmax %.10f\n",ortmax);
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
      printf("debug: after optimaztion ortmax %.10f\n",ortmax);
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
//using a dot product, and then determine whether they align with the expected direction specified for SOL, PFR, or CORE.
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

  double p3[2]={gridzone->start_point_R[start], gridzone->first_bnd_curve->points[start].y};
  double p4[2]={gridzone->start_point_R[start+1], gridzone->first_bnd_curve->points[start+1].y};

  if(dot_product(p2[0]-p1[0], p2[1]-p1[1], p4[0]-p3[0], p4[1]-p3[1])*dir<0.0)
  {
    func->rescale[0]=-1.0;
    func->rescale[1]=-1.0;
    printf("DEBUG reverse the poloidal diretction in %s.\n",gridzone->name);
  }
}

// tracing from a grid zone start point and finally arrive to the end curve.
// Here we assume the tracing direction HAS BEEN CORRECTED by change_poloidal_direction
// We also asumme that the tracing line will arrive at the end curve.
// after use, do not forget to FREE it!!!!!
// We also asumme there is only ONE intersection on the end curve.
static Curve* create_GridTubeCurve_by_tracing(double start_R, double start_Z, Curve* end_curve, 
                                              ode_function* func,ode_solver* solver)
{
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
      // new_c->n_point>2 otherwise for the core region the frist element of curve is always
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
    if(new_c->n_point>50000)
    {
      printf("The points of the curve has %zu\n", new_c->n_point);
      fprintf(stderr,"WARNING: please check the points number in the curve.\n");
      exit(EXIT_FAILURE);
    }
  }
  // write_curve("DEBUG_tracing_c", new_c);
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

  Curve* curr_c=create_GridTubeCurve_by_tracing(gridzone->start_point_R[1],
                                                gridzone->start_point_Z[1],
                                                gridzone->end_curve,
                                                func, solver);
  Curve* curr_gpt_c=create_curve(np);

  for(int i=0;i<np;i++)
  {
    add_last_point_curve(curr_gpt_c, 0.0, 0.0);
  }
  double* len_curr_gpt_c=malloc(np*sizeof(double));
  GirdTubeStr* gridtube=create_GridTube(prev_c, prev_gpt_c, len_prev_gpt_c, 
                                        curr_c, curr_gpt_c, len_curr_gpt_c,
                                        gridzone->guard_start[0],
                                        gridzone->guard_end[0],
                                        gridzone->pasmin[0]);

  // printf("%s %s %s\n","guard_start", "guard_end","pasmin");
  // printf("%lf %lf %lf\n", gridzone->guard_start[0], gridzone->guard_end[0],gridzone->pasmin[0]);

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

    curr_c=create_GridTubeCurve_by_tracing(gridzone->start_point_R[i],
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

    calc_points_from_CARRE(gridtube);
  
    sprintf(name,"%s%d",gridzone->name, i);
    write_curve(name,curr_gpt_c);
    
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


