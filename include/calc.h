#ifndef CALC_h
#define CALC_h



int VIntersect(double x1,double y1,double x2,double y2,
    double x3,double y3,double x4,double y4,double* ar,double* br);
#endif


// Macros from DivGeo
#ifndef min
#define min(x,y) ((x)<(y) ? (x) : (y))
#endif
#ifndef max
#define max(x,y) ((x)>(y) ? (x) : (y))
#endif
#ifndef inrange_s
#define inrange_s(a,b1,b2) ((a)<max((b1),(b2)) && (a)>min((b1),(b2)))
#endif
#ifndef swap
#define swap(x,y) ((x)+=(y),(y)=(x)-(y),(x)-=(y))
#endif
