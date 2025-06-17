#include "twodimgridgen.h"
#include "carrefunction.h"
#include <math.h>
#define DEFAULT_MARGIN 20

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

void generate_CARRE_2Dgrid(TwoDimGrid* grid,
                           GridZone* gridzone,
                           ode_function* func,
                           ode_solver* solver)
{
  return;
}