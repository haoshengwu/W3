#include <stdlib.h>
#include <math.h>
#include "calc.h"


// refer to DivGeo [calc.c]
/* Return -1 on non-intersection, 0 on intersection */
/* (x1,y1,x2,y2) is the first segment, and (x3,y3,x4,y4) is the second.
   *ar and *br are how far along the intersection is on each segment. */

int VIntersect(double x1,double y1,double x2,double y2, \
double x3,double y3,double x4,double y4,double* ar,double* br) 
{
  double a,b,d;
  d=(y4-y3)*(x2-x1)-(x4-x3)*(y2-y1);
  if (d==0) return -2;
  a=(y4-y3)*(x3-x1)-(x4-x3)*(y3-y1);
  b=(y2-y1)*(x3-x1)-(x2-x1)*(y3-y1);
  a=a/d;
  b=b/d;

  if (ar!=NULL) *ar=a;
  if (br!=NULL) *br=b;
  if (a<0 || a>1 || b<0 || b>1) return -1;
  return 0;
}

void bilinear_coeff(double** f, double* x, double* y,
                    const int nx, const int ny,
                    double*** coeff)
{
    // 遍历每个网格单元，计算插值系数
    for (int i = 0; i < nx - 1; i++) {
        for (int j = 0; j < ny - 1; j++) {
            double dx = x[i + 1] - x[i];
            double dy = y[j + 1] - y[j];
            double deltaxy = dx * dy;

            // 简化中间变量
            double f00 = f[i][j], f10 = f[i+1][j];
            double f01 = f[i][j+1], f11 = f[i+1][j+1];

            // 计算插值系数
            coeff[i][j][0] = (f00 * x[i+1] * y[j+1] - f10 * x[i] * y[j+1]
                            - f01 * x[i+1] * y[j] + f11 * x[i] * y[j]) / deltaxy;

            coeff[i][j][1] = (-f00 * y[j+1] + f10 * y[j+1]
                             + f01 * y[j] - f11 * y[j]) / deltaxy;

            coeff[i][j][2] = (-f00 * x[i+1] + f10 * x[i]
                             + f01 * x[i+1] - f11 * x[i]) / deltaxy;

            coeff[i][j][3] = (f00 - f10 - f01 + f11) / deltaxy;
    }
}
}

void rphi_to_XY(const double r, const double phi, double *x, double *y)
{
  double phi_rad = (phi) * PI / 180.0;
  *x = r * cos(phi_rad);
  *y = r * sin(phi_rad);
}

int ifind(double x, double* t, int n, int inc) 
{
  // If x is less than or equal to the first element, return the first index
  if (x <= t[0]) 
  {
    return 0; // Index in C is 0-based
  }
  // If x is greater than or equal to the last segment, return the last valid index
  else if (x >= t[n - inc]) 
  {
    return n - 2 * inc; // Adjusted for 0-based indexing
  } 
  else 
  {
    int ifind = 0;
    int j = n / inc; // Calculate the number of segments
    // Binary search for the correct index
    while (j - ifind > 1) 
    {
      int jp = (j + ifind) / 2; // Midpoint
      if (x > t[inc * (jp - 1)]) 
      {
        ifind = jp; // Narrow the search to the upper half
      } 
      else 
      {
        j = jp; // Narrow the search to the lower half
      }
    }
    return inc * (ifind - 1); // Adjust for the step size and 0-based indexing
  }
}


void derivation(const double *f, const double *x, const int n, double *df)
/* f: the function f(x) value
   x: the x value
   n: the number of x
  df: the derivation of f(x), df(x)
*/
{
  for (int i=1; i<n-1;i++)
    {
      df[i] = ((f[i+1]-f[i])*(x[i]-x[i-1])/(x[i+1]-x[i])	\
	+(f[i]-f[i-1])*(x[i+1]-x[i])/(x[i]-x[i-1]))             \
       /(x[i+1]-x[i-1]);
     }
  double x2 = x[2]-x[1];
  double x3 = x[3]-x[1];
  double f2 = f[2]-f[1];
  double f3 = f[3]-f[1];
  df[0] = (f2*x3*x3-f3*x2*x2)/(x2*x3*(x3-x2));
  
  x2 = x[n-2]-x[n-1];
  x3 = x[n-3]-x[n-1];
  f2 = f[n-2]-f[n-1];
  f3 = f[n-3]-f[n-1];
 df[n-1] = (f2*x3*x3-f3*x2*x2)/(x2*x3*(x3-x2));
  
}

void derivation_2d(double **f, const double *x, const int nx,
                   const double *y, const int ny, double ***df)
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

  
