#include "mathbase.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define MAX_ITER_NEWTON 1000
#define epsilon 1.0E-12


void swap_double(double* a, double* b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}
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
        df[i][j][1] = (-f[i][j+2] + 8 * f[i][j+1] - 8 * f[i][j-1] + f[i][j-2])/( 12 * dy);
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
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy)
{
  // assume uniform dx and dy
  double dx = (x[nx-1] - x[0])/(nx-1);
  double dy = (y[ny-1] - y[0])/(ny-1);
  double dxdy = dx*dy;
  //precision cretirea
  double eps = fmin(dx, dy) * 1e-6;
  
  //examine the boundary
  if (target_x < x[0] - eps || target_x > x[nx-1] + eps) 
  {
    fprintf(stderr, "Error: target_x (%.6f) is out of bounds [%.12f, %.12f]\n", target_x, x[0], x[nx-1]);
    exit(EXIT_FAILURE);
  }
  if (target_y < y[0] - eps || target_y > y[ny-1] + eps) 
  {
    fprintf(stderr, "Error: target_y (%.6f) is out of bounds [%.12f, %.12f]\n", target_y, y[0], y[ny-1]);
    exit(EXIT_FAILURE);
  }


  if (fabs(target_x - x[0]) < eps) target_x = x[0];
  if (fabs(target_x - x[nx-1]) < eps) target_x = x[nx-1];
  if (fabs(target_y - y[0]) < eps) target_y = y[0];
  if (fabs(target_y - y[ny-1]) < eps) target_y = y[ny-1];

  //nx points then nx-1 cells, ny points then ny-1 cells
  int i = floor((target_x - x[0])/dx);
  int j = floor((target_y - y[0])/dy);

  if (i < 0) i = 0;
  if (i >= nx - 1) i = nx - 2;
  if (j < 0) j = 0;
  if (j >= ny - 1) j = ny - 2;
  
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
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy)
{
  if (dfdx == NULL || dfdy == NULL || d2fdxdy == NULL)
  {
    printf("df/dx, df/dy, d2f/dxdy are calculted by f\n");
  }
  else 
  {
    printf("df/dx, df/dy, d2f/dxdy are provided by user\n");
  }
  return;
}

void cubicherm_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy)
{
  //This function is according to pspline function dnherm2() and function herm2fcn()
  
  // assume uniform dx and dy
  double dx = (x[nx-1] - x[0])/(nx-1);
  double dy = (y[ny-1] - y[0])/(ny-1);
  double dxdy = dx*dy;
  //precision cretirea
  double eps = fmin(dx, dy) * 1e-6;
  
  //examine the boundary
  if (target_x < x[0] - eps || target_x > x[nx-1] + eps) 
  {
    fprintf(stderr, "Error: target_x (%.6f) is out of bounds [%.12f, %.12f]\n", target_x, x[0], x[nx-1]);
    exit(EXIT_FAILURE);
  }
  if (target_y < y[0] - eps || target_y > y[ny-1] + eps) 
  {
    fprintf(stderr, "Error: target_y (%.6f) is out of bounds [%.12f, %.12f]\n", target_y, y[0], y[ny-1]);
    exit(EXIT_FAILURE);
  }

  if (fabs(target_x - x[0]) < eps) target_x = x[0];
  if (fabs(target_x - x[nx-1]) < eps) target_x = x[nx-1];
  if (fabs(target_y - y[0]) < eps) target_y = y[0];
  if (fabs(target_y - y[ny-1]) < eps) target_y = y[ny-1];

  //nx points then nx-1 cells, ny points then ny-1 cells
  int xc = floor((target_x - x[0])/dx);
  int yc = floor((target_y - y[0])/dy);

  if (xc < 0) xc = 0;
  if (xc >= nx - 1) xc = nx - 2;
  if (yc < 0) yc = 0;
  if (yc >= ny - 1) yc = ny - 2;


  double fxy_tmp[2][2][2];
  double dfdx_tmp[2][2][2];
  double dfdy_tmp[2][2][2];
  double d2fdxdy_tmp[2][2][2];

  if (dfdx == NULL || dfdy == NULL || d2fdxdy == NULL )
  {
    //printf("df/dx, df/dy, d2f/dxdy are calculted by f\n");
    for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
      {
        int iyp=min(ny-1,yc+j+1);
        int iym=max(0,yc+j-1);
        int ixp=min(nx-1,xc+i+1);
        int ixm=max(0,xc+i-1);

        for (int k=0; k<2; k++)
        {
        fxy_tmp[i][j][k] = f[xc+i][yc+j][k];
        dfdx_tmp[i][j][k] = (f[ixp][yc+j][k] - f[ixm][yc+j][k])/(x[ixp]-x[ixm]);
        dfdy_tmp[i][j][k] = (f[xc+i][iyp][k] - f[xc+i][iym][k])/(y[iyp]-y[iym]);
        d2fdxdy_tmp[i][j][k] = (f[ixp][iyp][k] - f[ixm][iyp][k] - f[ixp][iym][k] + f[ixm][iym][k])
                           / ((x[ixp]-x[ixm])*(y[iyp]-y[iym]));
        }
      }
    }
  }
  else 
  {
    for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
      {
        for (int k=0; k<2; k++)
        {
        fxy_tmp[i][j][k] = f[xc+i][yc+j][k];
        dfdx_tmp[i][j][k] = dfdx[xc+i][yc+j][k];
        dfdy_tmp[i][j][k] = dfdy[xc+i][yc+j][k];
        d2fdxdy_tmp[i][j][k] = d2fdxdy[xc+i][yc+j][k];
        }
      }
    }
  }
  
  double hx = dx;
  double hy = dy;
  double hxi = 1/hx;
  double hyi = 1/hy;

  double xp = (target_x-x[xc])*hxi;
  double xpi = 1-xp;
  double xp2 = xp*xp;
  double xpi2 = xpi*xpi;
  double ax = xp2*(3.0-2.0*xp);
  double axbar = 1.0-ax;
  double bx = -xp2*xpi;
  double bxbar = xpi2*xp;

  double yp = (target_y-y[yc])*hyi;
  double ypi = 1-yp;
  double yp2 = yp*yp;
  double ypi2 = ypi*ypi;
  double ay = yp2*(3.0-2.0*yp);
  double aybar = 1.0-ay;
  double by = -yp2*ypi;
  double bybar = ypi2*yp;

  double axp=6.0*xp*xpi;
  double axbarp=-axp;
  double bxp=xp*(3.0*xp-2.0);
  double bxbarp=xpi*(3.0*xpi-2.0);

  double ayp=6.0*yp*ypi;
  double aybarp=-ayp;
  double byp=yp*(3.0*yp-2.0);
  double bybarp=ypi*(3.0*ypi-2.0);

  double sum[2]={0.0, 0.0};

  //printf("xp: %lf\n",xp);
  //printf("yp: %lf\n",yp);
  for (int k=0;k<2;k++)
  {

    sum[k] = axbar * (aybar * fxy_tmp[0][0][k] + ay * fxy_tmp[0][1][k]) +
                ax * (aybar * fxy_tmp[1][0][k] + ay * fxy_tmp[1][1][k]);

    sum[k] += hx * (bxbar * (aybar * dfdx_tmp[0][0][k] + ay * dfdx_tmp[0][1][k]) +
                       bx * (aybar * dfdx_tmp[1][0][k] + ay * dfdx_tmp[1][1][k]));

    sum[k] += hy * (axbar * (bybar * dfdy_tmp[0][0][k] + by * dfdy_tmp[0][1][k]) +
                       ax * (bybar * dfdy_tmp[1][0][k] + by * dfdy_tmp[1][1][k]));

    sum[k] += hx * hy * (bxbar * (bybar * d2fdxdy_tmp[0][0][k] + by * d2fdxdy_tmp[0][1][k]) +
                            bx * (bybar * d2fdxdy_tmp[1][0][k] + by * d2fdxdy_tmp[1][1][k]));
  }
  *value1 = sum[0];
  *value2 = sum[1];
  return;
}



void bilenar_1d(double target_x, double target_y, int nx, double *x,  int ny, double *y, 
                double **f, double *value, double **dfdx, double **dfdy, double **d2fdxdy)
{
// assume uniform dx and dy
  double dx = (x[nx-1] - x[0])/(nx-1);
  double dy = (y[ny-1] - y[0])/(ny-1);
  double dxdy = dx*dy;
  //precision cretirea
  double eps = fmin(dx, dy) * 1e-6;
  
  //examine the boundary
  if (target_x < x[0] - eps || target_x > x[nx-1] + eps) 
  {
    fprintf(stderr, "Error: target_x (%.6f) is out of bounds [%.12f, %.12f]\n", target_x, x[0], x[nx-1]);
    exit(EXIT_FAILURE);
  }
  if (target_y < y[0] - eps || target_y > y[ny-1] + eps) 
  {
    fprintf(stderr, "Error: target_y (%.6f) is out of bounds [%.12f, %.12f]\n", target_y, y[0], y[ny-1]);
    exit(EXIT_FAILURE);
  }


  if (fabs(target_x - x[0]) < eps) target_x = x[0];
  if (fabs(target_x - x[nx-1]) < eps) target_x = x[nx-1];
  if (fabs(target_y - y[0]) < eps) target_y = y[0];
  if (fabs(target_y - y[ny-1]) < eps) target_y = y[ny-1];

  //nx points then nx-1 cells, ny points then ny-1 cells
  int i = floor((target_x - x[0])/dx);
  int j = floor((target_y - y[0])/dy);

  if (i < 0) i = 0;
  if (i >= nx - 1) i = nx - 2;
  if (j < 0) j = 0;
  if (j >= ny - 1) j = ny - 2;
  
  double f00 = f[i][j];
  double f10 = f[i+1][j];
  double f01 = f[i][j+1];
  double f11 = f[i+1][j+1];


  double x0 = x[i], x1 = x[i+1];
  double y0 = y[j], y1 = y[j+1];
  double tx = (target_x - x0) / (x1 - x0);
  double ty = (target_y - y0) / (y1 - y0);
  
  *value = (1 - tx) * (1 - ty) * f00 +
            tx * (1 - ty) * f10 +
            (1 - tx) * ty * f01 +
            tx * ty * f11;
  return;           
}

void cubicherm_1d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double **f, double *value, double **dfdx, double **dfdy, double **d2fdxdy)
{
  //This function is according to pspline function dnherm2() and function herm2fcn()
  // assume uniform dx and dy
  double dx = (x[nx-1] - x[0])/(nx-1);
  double dy = (y[ny-1] - y[0])/(ny-1);
  double dxdy = dx*dy;
  //precision cretirea
  double eps = fmin(dx, dy) * 1e-6;
  
  //examine the boundary
  if (target_x < x[0] - eps || target_x > x[nx-1] + eps) 
  {
    fprintf(stderr, "Error: target_x (%.6f) is out of bounds [%.12f, %.12f]\n", target_x, x[0], x[nx-1]);
    exit(EXIT_FAILURE);
  }
  if (target_y < y[0] - eps || target_y > y[ny-1] + eps) 
  {
    fprintf(stderr, "Error: target_y (%.6f) is out of bounds [%.12f, %.12f]\n", target_y, y[0], y[ny-1]);
    exit(EXIT_FAILURE);
  }

  if (fabs(target_x - x[0]) < eps) target_x = x[0];
  if (fabs(target_x - x[nx-1]) < eps) target_x = x[nx-1];
  if (fabs(target_y - y[0]) < eps) target_y = y[0];
  if (fabs(target_y - y[ny-1]) < eps) target_y = y[ny-1];

  //nx points then nx-1 cells, ny points then ny-1 cells
  int xc = floor((target_x - x[0])/dx);
  int yc = floor((target_y - y[0])/dy);
  
  //nx points then nx-1 cells, ny points then ny-1 cells
  if (xc < 0) xc = 0;
  if (xc >= nx - 1) xc = nx - 2;
  if (yc < 0) yc = 0;
  if (yc >= ny - 1) yc = ny - 2;

  double fxy_tmp[2][2];
  double dfdx_tmp[2][2];
  double dfdy_tmp[2][2];
  double d2fdxdy_tmp[2][2];

  if (dfdx == NULL || dfdy == NULL || d2fdxdy == NULL )
  {
    //printf("df/dx, df/dy, d2f/dxdy are calculted by f\n");
    for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
      {
        int iyp=min(ny-1,yc+j+1);
        int iym=max(0,yc+j-1);
        int ixp=min(nx-1,xc+i+1);
        int ixm=max(0,xc+i-1);


        fxy_tmp[i][j] = f[xc+i][yc+j];
        dfdx_tmp[i][j] = (f[ixp][yc+j] - f[ixm][yc+j])/(x[ixp]-x[ixm]);
        dfdy_tmp[i][j]= (f[xc+i][iyp] - f[xc+i][iym])/(y[iyp]-y[iym]);
        d2fdxdy_tmp[i][j] = (f[ixp][iyp] - f[ixm][iyp] - f[ixp][iym] + f[ixm][iym])
                           / ((x[ixp]-x[ixm])*(y[iyp]-y[iym]));
      }
    }
  }
  else 
  {
    for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
      {
        fxy_tmp[i][j] = f[xc+i][yc+j];
        dfdx_tmp[i][j] = dfdx[xc+i][yc+j];
        dfdy_tmp[i][j] = dfdy[xc+i][yc+j];
        d2fdxdy_tmp[i][j]= d2fdxdy[xc+i][yc+j];
      }
    }
  }
  
  double hx = dx;
  double hy = dy;
  double hxi = 1.0/hx;
  double hyi = 1.0/hy;

  double xp = (target_x-x[xc])*hxi;
  double xpi = 1-xp;
  double xp2 = xp*xp;
  double xpi2 = xpi*xpi;
  double ax = xp2*(3.0-2.0*xp);
  double axbar = 1.0-ax;
  double bx = -xp2*xpi;
  double bxbar = xpi2*xp;

  double yp = (target_y-y[yc])*hyi;
  double ypi = 1-yp;
  double yp2 = yp*yp;
  double ypi2 = ypi*ypi;
  double ay = yp2*(3.0-2.0*yp);
  double aybar = 1.0-ay;
  double by = -yp2*ypi;
  double bybar = ypi2*yp;

  double axp=6.0*xp*xpi;
  double axbarp=-axp;
  double bxp=xp*(3.0*xp-2.0);
  double bxbarp=xpi*(3.0*xpi-2.0);

  double ayp=6.0*yp*ypi;
  double aybarp=-ayp;
  double byp=yp*(3.0*yp-2.0);
  double bybarp=ypi*(3.0*ypi-2.0);

  double sum=0.0;

  //printf("xp: %lf\n",xp);
  //printf("yp: %lf\n",yp);

  sum = axbar * (aybar * fxy_tmp[0][0] + ay * fxy_tmp[0][1]) +
              ax * (aybar * fxy_tmp[1][0]+ ay * fxy_tmp[1][1]);

  sum += hx * (bxbar * (aybar * dfdx_tmp[0][0] + ay * dfdx_tmp[0][1]) +
                     bx * (aybar * dfdx_tmp[1][0]+ ay * dfdx_tmp[1][1]));

  sum += hy * (axbar * (bybar * dfdy_tmp[0][0] + by * dfdy_tmp[0][1]) +
                     ax * (bybar * dfdy_tmp[1][0] + by * dfdy_tmp[1][1]));

  sum += hx * hy * (bxbar * (bybar * d2fdxdy_tmp[0][0] + by * d2fdxdy_tmp[0][1]) +
                          bx * (bybar * d2fdxdy_tmp[1][0] + by * d2fdxdy_tmp[1][1]));

  *value = sum;
  return;
}

// referen to herm1ev.f90 in pspline
void cubicherm1D_eval(const void* interp_data, double target_x, double* value){
  const CubicHerm1dData* data = (const CubicHerm1dData*) interp_data;
  int nx = data->nx;
// range check
  double dx = (data->x[nx-1] - data->x[0])/(nx-1);
  //precision cretirea
  double eps = dx * 1e-6;

  //examine the boundary
  if (target_x < data->x[0] - eps || target_x > data->x[nx-1] + eps) 
  {
    fprintf(stderr, "Error: target_x (%.6f) is out of bounds [%.12f, %.12f]\n", target_x, data->x[0], data->x[nx-1]);
    exit(EXIT_FAILURE);
  }
  if (fabs(target_x - data->x[0]) < eps) target_x = data->x[0];
  if (fabs(target_x - data->x[nx-1]) < eps) target_x = data->x[nx-1];
  //nx points then nx-1 cells, ny points then ny-1 cells
  int xc = floor((target_x - data->x[0])/dx);
  if (xc < 0) xc = 0;
  if (xc >= nx - 1) xc = nx - 2;

  double hx=(data->x[xc+1]-data->x[xc]);
  double hxi = 1.0/hx;

  double xp = (target_x-data->x[xc])*hxi;
  double xpi = 1-xp;
  double xp2 = xp*xp;
  double xpi2 = xpi*xpi;
  double ax = xp2*(3.0-2.0*xp);
  double axbar = 1.0-ax;
  double bx = -xp2*xpi;
  double bxbar = xpi2*xp;

  double sum=0.0;
  sum = axbar*data->fx[xc] + ax*data->fx[xc+1];
  sum = sum + hx*(bxbar*data->dfdx[xc]+bx*data->dfdx[xc+1]);

  *value = sum;
  return;
}


void cubicherm1D_deriv(const void* interp_data, double target_x, double* value)
{
  const CubicHerm1dData* data = (const CubicHerm1dData*) interp_data;
  int nx = data->nx;
// range check
  double dx = (data->x[nx-1] - data->x[0])/(nx-1);
  //precision cretirea
  double eps = dx * 1e-6;

  //examine the boundary
  if (target_x < data->x[0] - eps || target_x > data->x[nx-1] + eps) 
  {
    fprintf(stderr, "Error: target_x (%.6f) is out of bounds [%.12f, %.12f]\n", target_x, data->x[0], data->x[nx-1]);
    exit(EXIT_FAILURE);
  }
  if (fabs(target_x - data->x[0]) < eps) target_x = data->x[0];
  if (fabs(target_x - data->x[nx-1]) < eps) target_x = data->x[nx-1];
  //nx points then nx-1 cells, ny points then ny-1 cells
  int xc = floor((target_x - data->x[0])/dx);
  if (xc < 0) xc = 0;
  if (xc >= nx - 1) xc = nx - 2;
  double hx=(data->x[xc+1]-data->x[xc]);
  double hxi = 1.0/hx;

  double xp = (target_x-data->x[xc])*hxi;
  double xpi = 1-xp;
  double xp2 = xp*xp;
  double xpi2 = xpi*xpi;
  double axp=6.0*xp*xpi;
  double axbarp=-axp;
  double bxp=xp*(3.0*xp-2.0);
  double bxbarp=xpi*(3.0*xpi-2.0);
  double sum=0.0;
  sum=hxi*(axbarp*data->fx[xc] +axp*data->fx[xc+1]);
  sum=sum+bxbarp*data->fx[xc] + bxp*data->fx[xc+1];
  *value = sum;
  return;
}

Interp1DFunction* create_cubicherm1D_interp(double* x, double* fx, double* dfdx, int nx)
{
  CubicHerm1dData* data = malloc(sizeof(CubicHerm1dData));
  
  if (!data) {
    fprintf(stderr, "Failed to allocate CubicHerm1dData.\n");
    exit(EXIT_FAILURE);
  }

  bool input_xfxdfdx = false;

  data->nx = nx;
  data->x = malloc(sizeof(double) * nx);
  data->fx = malloc(sizeof(double) * nx);
  data->dfdx = malloc(sizeof(double) * nx);

  if (!data->x || !data->fx || !data->dfdx) {
    fprintf(stderr, "Failed to allocate arrays.\n");
    exit(EXIT_FAILURE);
  }

  if(x==NULL && fx==NULL && dfdx==NULL)
  {
    for(int i=0;i<nx;i++)
    {
      data->x[i]=0;
      data->fx[i]=0;
      data->dfdx[i]=0;
    }
  }
  else if(x!=NULL && fx!=NULL && dfdx!=NULL)
  {
    memcpy(data->x, x, sizeof(double) * nx);
    memcpy(data->fx, fx, sizeof(double) * nx);
    memcpy(data->dfdx, dfdx, sizeof(double) * nx);
  }
  else
  {
    free(data->x);
    free(data->fx);
    free(data->dfdx);
    fprintf(stderr, "Please check the input for the interpolator\n");
    exit(EXIT_FAILURE);
  }

  Interp1DFunction* interp = malloc(sizeof(Interp1DFunction));
  interp->data = data;
  interp->eval = cubicherm1D_eval;
  interp->deriv = cubicherm1D_deriv;
  interp->free_data = cubicherm1D_free_data;
  interp->name = "CubicHerm1D";
  return interp;
}

//modify the CubicHerm1dData value! not the addreass. THE SIZE should be exactly the same.
void modify_cubicherm1D_data(CubicHerm1dData* data, double* x, double* fx, double* dfdx, int nx)
{
  if(nx!=data->nx)
  {
    printf("The size of new data is not same with previous!!!\n");
    exit(EXIT_FAILURE);
  }
  // printf("DEBUG in modify_cubicherm1D_data\n");
  for(int i=0;i<nx;i++)
  {
    if(x!=NULL) data->x[i]=x[i];
    if(fx!=NULL) data->fx[i]=fx[i];
    if(dfdx!=NULL)data->dfdx[i]=dfdx[i];
  }
  return;
}

void cubicherm1D_free_data(void** ptr)
{
  if (ptr && *ptr)
  {
    CubicHerm1dData* data = (CubicHerm1dData*)(*ptr);
      free(data->x);
      free(data->fx);
      free(data->dfdx);
      free(data); 
      *ptr = NULL;
  }
}
void free_interp1d_function(Interp1DFunction* interp)
{
  if (interp == NULL) return;
  if (interp->free_data && interp->data) 
    {
      interp->free_data(&(interp->data));
    }
    free(interp);
    printf("Finish free interpolator.\n");
}



//assume cubichermit interpolation, will be update in future.
//Use newton method to calculate x which have the fx value;
void Newton_Raphson_Method(double* x, const double* fx_target, Interp1DFunction* interp)
{
  double x_prev, x_curr;
  CubicHerm1dData* data = (CubicHerm1dData*) interp->data;
  int nx = data->nx;
  x_prev=data->x[0]+(data->x[nx-1]-data->x[0])*1.0e-6;
  double fx_prev, dfdx_prev;
  double fx_curr;
  for(int i=0; i<MAX_ITER_NEWTON; i++)
  {
    interp->eval(interp->data, x_prev, &fx_prev);
    interp->deriv(interp->data, x_prev, &dfdx_prev);
    if (fabs(dfdx_prev) < 1e-14) 
    {
      printf("Warning: derivative too small in Newton iteration.\n");
      break;
    }
    x_curr=x_prev-(fx_prev-*fx_target)/dfdx_prev;
    interp->eval(interp->data, x_curr, &fx_curr);
    if(fabs(fx_curr-*fx_target)<epsilon)
    {
      *x = x_curr;
      return;
    }
    x_prev = x_curr;
  }
  printf("Newton_Raphson_Method: maximum iterations reached.\n");
  *x = x_curr;
  return;
}

//assume cubichermit interpolation, will be update in future.
void Secant_Method(double* x, const double* fx_target, Interp1DFunction* interp)
{
  double x_prev, x_curr, x_next;
  int nx;
  if(strcmp(interp->name, "CubicHerm1D")==0)
  {
    CubicHerm1dData* data = (CubicHerm1dData*) interp->data;
    nx = data->nx;
    x_prev = data->x[0] ;
    x_curr = data->x[nx - 1];
  }
  else
  {
    printf("Only support CubicHermit data!\n");
    exit(EXIT_FAILURE);
  }

  double fx_prev, fx_curr, fx_next;

  interp->eval(interp->data, x_prev, &fx_prev);
  interp->eval(interp->data, x_curr, &fx_curr);

  for (int i = 0; i < MAX_ITER_NEWTON; ++i) 
  {
    double denom = fx_curr - fx_prev;
    if (fabs(denom) < 1e-14) {
      printf("Warning: function difference too small in Secant iteration.\n");
      break;
    }

    double x_next = x_curr - (fx_curr - *fx_target) * (x_curr - x_prev) / denom;

    interp->eval(interp->data, x_next, &fx_next);

    if (fabs(fx_next - *fx_target) < epsilon) {
      *x = x_next;
      return;
    }

    x_prev = x_curr;
    fx_prev = fx_curr;
    x_curr = x_next;
    fx_curr = fx_next;
  }

  printf("Secant_Method: maximum iterations reached.\n");
  *x = x_curr;  // fallback
}