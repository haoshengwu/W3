#include <stdlib.h>

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

void bilinear_coeff(const double** f, const double *x, const double *y,
                   const int* nx, const int* ny,
                   double*** coeff)
{

  for (int i=0; i < *nx-1; i++)
    for (int j = 0; j < *ny-1; j++)
    {
      double deltaxy = (x[i+1]-x[i])*(y[j+1]-y[j]);
      coeff[i][j][0] = (f[i][j]*x[i]*y[i] - f[i+1][j]*x[i]*y[j+1] \
        - f[i][j+1]*x[i+1]*y[j] + f[i+1][j+1]*x[i]*y[j])/deltaxy;
      
      coeff[i][j][1] = (-f[i][j]*y[j+1] + f[i+1][j]*y[j+1] \
                        + f[i][j+1]*y[j] - f[i+1][j+1]*y[j])/deltaxy;
      
      coeff[i][j][2] = (-f[i][j]*x[j+1] + f[i+1][j]*x[i]	\
			+ f[i][j+1]*y[i+1] - f[i+1][j+1]*x[i])/deltaxy;

      coeff[i][j][3] = (f[i][j] - f[i+1][j] - f[i][j+1] \ 
			+ f[i+1][j+1])/deltaxy;
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
      df[i][j][0] = ((f[i + 1][j] - f[i][j]) * (x[i] - x[i - 1]) / (x[i + 1] - x[i]) +
                    (f[i][j] - f[i - 1][j]) * (x[i + 1] - x[i]) / (x[i] - x[i - 1])) /
                    (x[i + 1] - x[i - 1]);

      df[i][j][1] = ((f[i][j + 1] - f[i][j]) * (y[j] - y[j - 1]) / (y[j + 1] - y[j]) +
                     (f[i][j] - f[i][j - 1]) * (y[j + 1] - y[j]) / (y[j] - y[j - 1])) /
                     (y[j + 1] - y[j - 1]);
    }
  }

  for (int j = 1; j < ny - 1; j++)
  {
    // i = 0 (left boundary))
    double x2 = x[1] - x[0];
    double x3 = x[2] - x[0];
    double f2 = f[1][j] - f[0][j];
    double f3 = f[2][j] - f[0][j];
    df[0][j][0] = (f2 * x3 * x3 - f3 * x2 * x2) / (x2 * x3 * (x3 - x2));

    // i = n-1 (right boundary)
    x2 = x[nx - 2] - x[nx - 1];
    x3 = x[nx - 3] - x[nx - 1];
    f2 = f[nx - 2][j] - f[nx - 1][j];
    f3 = f[nx - 3][j] - f[nx - 1][j];
    df[nx - 1][j][0] = (f2 * x3 * x3 - f3 * x2 * x2) / (x2 * x3 * (x3 - x2));
  }
  for (int i = 1; i < nx - 1; i++) 
  {
  // j = 0 (bottom boundary)
    double y2 = y[1] - y[0];
    double y3 = y[2] - y[0];
    double f2 = f[i][1] - f[i][0];
    double f3 = f[i][2] - f[i][0];
    df[i][0][1] = (f2 * y3 * y3 - f3 * y2 * y2) / (y2 * y3 * (y3 - y2));

  // j = m-1 (top boundary)
    y2 = y[ny - 2] - y[ny - 1];
    y3 = y[ny - 3] - y[ny - 1];
    f2 = f[i][ny - 2] - f[i][ny - 1];
    f3 = f[i][ny - 3] - f[i][ny - 1];
    df[i][ny - 1][1] = (f2 * y3 * y3 - f3 * y2 * y2) / (y2 * y3 * (y3 - y2));
  }

  df[0][0][0] = (f[1][0] - f[0][0]) / (x[1] - x[0]);
  df[0][0][1] = (f[0][1] - f[0][0]) / (y[1] - y[0]);

  df[0][ny - 1][0] = (f[1][ny - 1] - f[0][ny - 1]) / (x[1] - x[0]);
  df[0][ny - 1][1] = (f[0][ny - 1] - f[0][ny - 2]) / (y[ny - 1] - y[ny - 2]);

  df[nx - 1][0][0] = (f[nx - 1][0] - f[nx - 2][0]) / (x[nx - 1] - x[nx - 2]);
  df[nx - 1][0][1] = (f[nx - 1][1] - f[nx - 1][0]) / (y[1] - y[0]);

  df[nx - 1][ny - 1][0] = (f[nx - 1][ny - 1] - f[nx - 2][ny - 1]) / (x[nx - 1] - x[nx - 2]);
  df[nx - 1][ny - 1][1] = (f[nx - 1][ny - 1] - f[nx - 1][ny - 2]) / (y[ny - 1] - y[ny - 2]);

}

  
