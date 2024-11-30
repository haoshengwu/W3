#include "mathbase.h"


double deg2rad(double phi)
{
  return phi * M_PI / 180.0;
}

void central_2nd_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df)
{
  for (int i = 1; i < nx - 1; i++) 
  {
    for (int j = 1; j < ny - 1; j++)
    {
      df[i][j][0] = (f[i+1][j] - f[i-1][j]) / (x[i+1] - x[i-1]);
      df[i][j][1] = (f[i][j+1] - f[i][j-1]) / (y[j+1] - y[j-1]);
    }
  }

  for (int j = 1; j < ny - 1; j++)
  {
    df[0][j][0] = (f[1][j] - f[0][j]) / (x[1] - x[0]);
    df[0][j][1] = (f[0][j+1] - f[0][j-1]) / (y[j+1] - y[j-1]);
  }

  for (int j = 1; j < ny - 1; j++) 
  {
    df[nx-1][j][0] = (f[nx-1][j] - f[nx-2][j]) / (x[nx-1] - x[nx-2]);
    df[nx-1][j][1] = (f[nx-1][j+1] - f[nx-1][j-1]) / (y[j+1] - y[j-1]);
  }

  for (int i = 1; i < nx - 1; i++)
  {
    df[i][0][0] = (f[i+1][0] - f[i-1][0]) / (x[i+1] - x[i-1]);
    df[i][0][1] = (f[i][1] - f[i][0]) / (y[1] - y[0]);
  }

  for (int i = 1; i < nx - 1; i++)
  {
    df[i][ny-1][0] = (f[i+1][ny-1] - f[i-1][ny-1]) / (x[i+1] - x[i-1]);
    df[i][ny-1][1] = (f[i][ny-1] - f[i][ny-2]) / (y[ny-1] - y[ny-2]);
  }

    df[0][0][0] = (f[1][0] - f[0][0]) / (x[1] - x[0]);
    df[0][0][1] = (f[0][1] - f[0][0]) / (y[1] - y[0]);

    df[0][ny-1][0] = (f[1][ny-1] - f[0][ny-1]) / (x[1] - x[0]);
    df[0][ny-1][1] = (f[0][ny-1] - f[0][ny-2]) / (y[ny-1] - y[ny-2]);

    df[nx-1][0][0] = (f[nx-1][0] - f[nx-2][0]) / (x[nx-1] - x[nx-2]);
    df[nx-1][0][1] = (f[nx-1][1] - f[nx-1][0]) / (y[1] - y[0]);

    df[nx-1][ny-1][0] = (f[nx-1][ny-1] - f[nx-2][ny-1]) / (x[nx-1] - x[nx-2]);
    df[nx-1][ny-1][1] = (f[nx-1][ny-1] - f[nx-1][ny-2]) / (y[ny-1] - y[ny-2]);
}  

void central_4th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df)
{
  if (nx < 5 || ny < 5) 
  {
    printf("ERROR: Grid size too small for 4th-order central difference!\n");
    return;
  }
  //assume uniform grid in x and y.
  double dx = x[1] - x[0];
  double dy = y[1] - y[0];
  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      if (i>=2 && i<=nx-3)
      {
        df[i][j][0] = (-f[i+2][j] + 8 * f[i+1][j] - 8 * f[i-1][j] + f[i-2][j])/(12 * dx);
      }
      else if (i<2)
      {
        df[i][j][0] = (-3 * f[i][j] + 4 * f[i+1][j] - f[i+2][j]) / (2 * dx);
      }
      else if (i>nx-3)
      {
        df[i][j][0] = (3 * f[i][j] - 4 * f[i-1][j] + f[i-2][j]) / (2 * dx);
      }

      if (j>=2 && j<=ny-3)
      {
        df[i][j][1] = (-f[i][j+2] + 8 * f[i][j+1] - 8 * f[i][j-1] + f[i][j-2])/(3 * (12 * dy));
      }
      else if (j<2)
      {
        df[i][j][1] = (-3 * f[i][j] + 4 * f[i][j+1] - f[i][j+2]) / (2 * dy);
      }
      else if (j>ny-3)
      {
        df[i][j][1] = (3 * f[i][j] - 4 * f[i][j-1] + f[i][j-2]) / (2 * dy);
      }
    }
  }
}

void bilenar_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y, 
                double ***f, double *value1, double *value2)
{
  // assume uniform dx and dy
  double dx = (x[nx-1] - x[0])/(nx-1);
  double dy = (y[ny-1] - y[0])/(ny-1);
  double dxdy = dx*dy;
  //nx points then nx-1 cells, ny points then ny-1 cells
  int i = floor((target_x - x[0])/dx);
  int j = floor((target_y - y[0])/dy);
  //check range
  if (i < 0 || i >= nx - 1 || j < 0 || j >= ny - 1) {
    printf("Target point is out of x\n");
    exit(1);
  }
  double f00 = f[i][j][0];
  double f10 = f[i+1][j][0];
  double f01 = f[i][j+1][0];
  double f11 = f[i+1][j+1][0];

  double g00 = f[i][j][1];
  double g10 = f[i+1][j][1];
  double g01 = f[i][j+1][1];
  double g11 = f[i+1][j+1][1];

  double x0 = x[i], x1 = x[i+1];
  double y0 = y[j], y1 = y[j+1];
  double tx = (target_x - x0) / (x1 - x0);
  double ty = (target_y - y0) / (y1 - y0);
  
  *value1 = (1 - tx) * (1 - ty) * f00 +
            tx * (1 - ty) * f10 +
            (1 - tx) * ty * f01 +
            tx * ty * f11;
  *value2 = (1 - tx) * (1 - ty) * g00 +
            tx * (1 - ty) * g10 +
            (1 - tx) * ty * g01 +
            tx * ty * g11;                  
}

void bicubic_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2)
{
  //TODO
  return;
}

void cubherm_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2)
{
  //TODO
  return;
}