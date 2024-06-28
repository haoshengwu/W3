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