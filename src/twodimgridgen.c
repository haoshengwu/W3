#include "twodimgridgen.h"
#include "carrefunction.h"
#include <math.h>
#define DEFAULT_MARGIN 20
#define MAX_NUM_TRACING 40000

typedef struct
{
/*********************************
* Input parameters
**********************************/
    //previous curve which is coresponding to (ir-1) gridcurve. 
    Curve *prev_c;

    //curve points which is coresponding to (ir )gridcurve.
    Curve *curr_c;
    
    //The length distributions of the grid pionts in the previous gridcurve. 
    double *len_prev_gridcurve;

    //the points' coordiantes for the grid which are along the prev_curve.
    Curve *prev_gridcurve;

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
    double *len_curr_gridcurve;
    Curve *curr_gridcurve;

} GirdTubeStr;

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
    double y2 = p2->y-base->y;
    double x2 = p2->x-base->x;
    double y1 = p1->y-base->y;
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

static void allocate_CarreOrthoProp(size_t n)
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
}

static void free_CarreOrthoProp(CarreOrthoStr *orthogonal)
{
  if(!orthogonal) return;
  free(orthogonal->ort);orthogonal->ort=NULL;
  free(orthogonal->ortpur); orthogonal->ortpur=NULL;
  free(orthogonal->propo); orthogonal->propo=NULL;
  free(orthogonal->varr); orthogonal->varr=NULL;
  free(orthogonal->tot); orthogonal->tot=NULL;
}

/*******************************************************************************
* This function calculate the orthogonalirty of the mesh points on two mesh curves.
* This fuction is refered to Fortran code CARRE/clort.F .
********************************************************************************/
//len_prev_gridcurve is length distribution for prev_gridcurve, 0.0 for the 1st and total length for the last
//Curve *prev_gridcurve is the previous curve contain the grid point
//len_curr_gridcurve is length distribution for curr_gridcurve, 0.0 for the 1st and total length for the last
//Curve *curr_gridcurve is the current curve contain the grid point
//'current curve' means we want the grid point in the current cure based on previous curve.
static void calc_ortho_from_CARRE(GirdTubeStr *gridtube, CarreOrthoStr *orthogonal)
{
  int n_point = gridtube->prev_gridcurve->n_point;
  if (n_point != gridtube->curr_gridcurve->n_point) {
    fprintf(stderr, "The points are not consistent in the two curves!\n");
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
  l1p  = calc_length(&gridtube->prev_gridcurve->points[0],
                     &gridtube->prev_gridcurve->points[1]);
  l2p  = calc_length(&gridtube->curr_gridcurve->points[0],
                     &gridtube->curr_gridcurve->points[1]);
  l12t = calc_length(&gridtube->prev_gridcurve->points[0],
                     &gridtube->curr_gridcurve->points[0]);

  double g1 = gridtube->guard_top;
  double g2 = gridtube->guard_end;

  for (int i = 1; i < n_point - 1; i++) 
  {
    double l12;
    if (g1 > 0.0 && g2 > 0.0) 
    {
      double fac1 = pow(g1 / (g1 + gridtube->len_prev_gridcurve[i]), 2.0);
      double fac2 = pow(g2 / (g2 + gridtube->len_prev_gridcurve[n_point - 1] -
                              gridtube->len_prev_gridcurve[i]), 2.0);
      double fac  = fac1 * (un - fac2) + fac2;
      l12 = fac * gridtube->len_prev_gridcurve[n_point - 1] +
            (1.0 - fac) * l12t;
    } 
    else 
    {
      l12 = l12t;
    }

    double l1m = l1p;
    double l2m = l2p;

    l1p = calc_length(&gridtube->prev_gridcurve->points[i],
                      &gridtube->prev_gridcurve->points[i + 1]);
    l2p = calc_length(&gridtube->curr_gridcurve->points[i],
                      &gridtube->curr_gridcurve->points[i + 1]);

    double cs1 = 0.0, cs2 = 0.0, cs3 = 0.0, cs4 = 0.0;
    if (i - 1 >= 0 && i + 1 < n_point) {
      cs1 = cosine_term(&gridtube->curr_gridcurve->points[i],
                        &gridtube->prev_gridcurve->points[i],
                        &gridtube->curr_gridcurve->points[i - 1]);
      cs2 = cosine_term(&gridtube->curr_gridcurve->points[i],
                        &gridtube->prev_gridcurve->points[i],
                        &gridtube->curr_gridcurve->points[i + 1]);
      cs3 = cosine_term(&gridtube->prev_gridcurve->points[i],
                        &gridtube->curr_gridcurve->points[i],
                        &gridtube->curr_gridcurve->points[i - 1]);
      cs4 = cosine_term(&gridtube->prev_gridcurve->points[i],
                        &gridtube->curr_gridcurve->points[i],
                        &gridtube->curr_gridcurve->points[i + 1]);
    } else 
    {
      printf("Warning: in calc_ortho_CARRE, Out of bounds access at i=%d\n", i);
    }

    cs1 /= (l2m * l12);
    cs2 /= (l2p * l12);
    cs3 /= (l1m * l12);
    cs4 /= (l1p * l12);

    double f1 = cs2 + cs3 - cs1 - cs4;
    double f2 = -(pow(g1 / gridtube->len_curr_gridcurve[i], 2.0) +
                  pow(g2 / (gridtube->len_curr_gridcurve[n_point - 1] -
                            gridtube->len_curr_gridcurve[i]), 2.0)) *
                (gridtube->len_curr_gridcurve[i] -
                 gridtube->len_prev_gridcurve[i] /
                 gridtube->len_prev_gridcurve[n_point - 1] *
                 gridtube->len_curr_gridcurve[n_point - 1]) /
                (gridtube->pasmin + g1 + g2);
    double f3 = pow(gridtube->pasmin /
                    (gridtube->len_curr_gridcurve[i] -
                     gridtube->len_curr_gridcurve[i - 1]), 2.0) -
                pow(gridtube->pasmin /
                    (gridtube->len_curr_gridcurve[i + 1] -
                     gridtube->len_curr_gridcurve[i]), 2.0);

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




static GirdTubeStr* allocate_GridTube(size_t n, double guard_top, double guard_end, double pasmin)
{
  if(guard_end<0.0||guard_end<0.0)
  {
    printf("WARNING: Guard length is less than ZERO!\n");
    printf("WARNING: Guard length is set to 0.0!\n");
    guard_end=0.0;
    guard_top=0.0;
  }
  GirdTubeStr* gridtube = malloc(sizeof(GirdTubeStr));
  if(!gridtube)
  {
    fprintf(stderr,"Failed to allocate gridtube.\n");
    exit(EXIT_FAILURE);
  }
  gridtube->pasmin=pasmin;
  gridtube->guard_top=guard_top;
  gridtube->guard_end=guard_end;
  gridtube->prev_c=create_curve(MAX_NUM_TRACING);
  gridtube->curr_c=create_curve(MAX_NUM_TRACING);
  gridtube->prev_gridcurve=create_curve(n);
  gridtube->curr_gridcurve=create_curve(n);
  gridtube->len_prev_gridcurve= calloc(n, sizeof(double));
  gridtube->len_curr_gridcurve= calloc(n, sizeof(double));
}

void free_GridTube(GirdTubeStr* gridtube)
{
  if (!gridtube) return;
  if (gridtube->prev_c) free_curve(gridtube->prev_c);
  if (gridtube->curr_c) free_curve(gridtube->curr_c);
  if (gridtube->prev_gridcurve) free_curve(gridtube->prev_gridcurve);
  if (gridtube->curr_gridcurve) free_curve(gridtube->curr_gridcurve);
  if (gridtube->len_prev_gridcurve) free(gridtube->len_prev_gridcurve);
  if (gridtube->len_curr_gridcurve) free(gridtube->len_curr_gridcurve);
  free(gridtube);
}



void generate_CARRE_2Dgrid_default(TwoDimGrid* grid,
                                   GridZoneInfo* gridzoneinfo,
                                   ode_function* func,
                                   ode_solver* solver)
{
  
  
  
  
  return;
}