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
                double ***f, double *value1, double *value2, void *intpl_data)
{
  // assume uniform dx and dy
  double dx = (x[nx-1] - x[0])/(nx-1);
  double dy = (y[ny-1] - y[0])/(ny-1);
  double dxdy = dx*dy;
  //nx points then nx-1 cells, ny points then ny-1 cells
  int i = floor((target_x - x[0])/dx);
  int j = floor((target_y - y[0])/dy);
  //check range
  //printf("debug: target_x: %lf, target_y: %lf\n", target_x, target_y);
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
                double ***f, double *value1, double *value2, void *intpl_data)
{
  if (intpl_data == NULL)
  {
    printf("df/dx, df/dy, d2f/dxdy are calculted by f\n");
  }
  else 
  {
    printf("df/dx, df/dy, d2f/dxdy are provided by user\n");
    bicubic_2d_data *pre_data = (bicubic_2d_data *) intpl_data;
  }
  return;
}

void cubicherm_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2, void *intpl_data)
{
  //This function is according to pspline function dnherm2() and function herm2fcn()
  // assume uniform dx and dy
  double dx = (x[nx-1] - x[0])/(nx-1);
  double dy = (y[ny-1] - y[0])/(ny-1);
  double dxdy = dx*dy;
  //nx points then nx-1 cells, ny points then ny-1 cells
  int xc = floor((target_x - x[0])/dx);
  int yc = floor((target_y - y[0])/dy);
  //check range
  //printf("debug: target_x: %lf, target_y: %lf\n", target_x, target_y);
  if (xc < 0 || xc >= nx - 1 || yc < 0 || yc >= ny - 1) 
  {
    printf("Target point is out of x\n");
    exit(1);
  }

  double fxy[2][2][2];
  double dfdx[2][2][2];
  double dfdy[2][2][2];
  double d2fdxdy[2][2][2];

  if (intpl_data == NULL)
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
        fxy[i][j][k] = f[xc+i][yc+j][k];
        dfdx[i][j][k] = (f[ixp][yc+j][k] - f[ixm][yc+j][k])/(x[ixp]-x[ixm]);
        dfdy[i][j][k] = (f[xc+i][iyp][k] - f[xc+i][iym][k])/(y[iyp]-y[iym]);
        d2fdxdy[i][j][k] = (f[ixp][iyp][k] - f[ixm][iyp][k] - f[ixp][iym][k] + f[ixm][iym][k])
                           / ((x[ixp]-x[ixm])*(y[iyp]-y[iym]));
        }
      }
    }
  }
  else 
  {
    printf("df/dx, df/dy, d2f/dxdy are provided by user\n");
    printf("it is not supported now\n");
    //todo
    exit(1);
    //bicubic_2d_data *pre_data = (bicubic_2d_data *) intpl_data;
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

    sum[k] = axbar * (aybar * fxy[0][0][k] + ay * fxy[0][1][k]) +
                ax * (aybar * fxy[1][0][k] + ay * fxy[1][1][k]);

    sum[k] += hx * (bxbar * (aybar * dfdx[0][0][k] + ay * dfdx[0][1][k]) +
                       bx * (aybar * dfdx[1][0][k] + ay * dfdx[1][1][k]));

    sum[k] += hy * (axbar * (bybar * dfdy[0][0][k] + by * dfdy[0][1][k]) +
                       ax * (bybar * dfdy[1][0][k] + by * dfdy[1][1][k]));

    sum[k] += hx * hy * (bxbar * (bybar * d2fdxdy[0][0][k] + by * d2fdxdy[0][1][k]) +
                            bx * (bybar * d2fdxdy[1][0][k] + by * d2fdxdy[1][1][k]));
  }
  *value1 = sum[0];
  *value2 = sum[1];
}