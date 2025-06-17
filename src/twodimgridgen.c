#include "twodimgridgen.h"
#include "carrefunction.h"

TwoDimGrid* create_2Dgrid(int npol, int nrad) 
{
  TwoDimGrid* grid = malloc(sizeof(TwoDimGrid));
  if (!grid) 
  {
        fprintf(stderr, "Failed to allocate TwoDimGrid.\n");
        exit(EXIT_FAILURE);
  }

  grid->npol = npol;
  grid->nrad = nrad;
  grid->points = malloc(sizeof(GridPoint) * npol * nrad);
  if (!grid->points) 
  {
        fprintf(stderr, "Failed to allocate grid points.\n");
        free(grid);
        exit(EXIT_FAILURE);
  }
  return grid;
}

void free_2Dgrid(TwoDimGrid* grid) 
{
  if (grid)
  {
    free(grid->points);
    free(grid);
  }
}

GridPoint get_point_2Dgrid(const TwoDimGrid* grid, int ip, int ir) 
{
  return grid->points[ir * grid->npol + ip];
}

double get_x_2Dgrid(const TwoDimGrid* grid, int ip, int ir) 
{
  return grid->points[ir * grid->npol + ip].x;
}

double get_y_2Dgrid(const TwoDimGrid* grid, int ip, int ir) 
{
  return grid->points[ir * grid->npol + ip].y;
}

void set_point_2Dgrid(TwoDimGrid* grid, int ir, int ip, double x, double y) 
{
  grid->points[ir * grid->npol + ip].x = x;
  grid->points[ir * grid->npol + ip].y = y;
}

void generate_CARRE_2Dgrid(TwoDimGrid* grid,
                           GridZone* gridzone,
                           ode_function* func,
                           ode_solver* solver)
{
  return;
}
