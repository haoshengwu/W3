#ifndef MATHBASE_H
#define MATHBASE_H

#define M_PI 3.14159265358979323846
#ifndef min
#define min(x,y) ((x)<(y) ? (x) : (y))
#endif
#ifndef max
#define max(x,y) ((x)>(y) ? (x) : (y))
#endif

void swap_double(double* a, double* b);

//convertion angle to radians
double deg2rad(double phi);

/******************************************************
 *    Difference Algorithm
 ******************************************************/
// df[ix][iy][0] is df/dx, is df[ix][iy][1] is df/dy
void central_2nd_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
void central_4th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
//ToDo
//void central_6th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);
//void central_8th_2d_diff(int nx, double *x,  int ny, double *y, double **f, double ***df);


/******************************************************
 *  Interpolation Algorithm
 *  'void *intpl_data' is used to point any intermediate structure data to speed up the calculation
 *******************************************************/
typedef struct 
{
  int nx;
  int ny;
  double **dfdx;
  double **dfdy;
  double **d2fdxdy;
} bicubic_2d_data;

typedef struct 
{
  double **dfdx;
  double **dfdy;
  double **d2fdxdy;
} CubicHerm2dData;

// dfdx dfdy d2fdxdy can be supported or just NULL pointer.
// here 2d means for TWO 2d-functions. NOT for ONE 2d-function.
void bilenar_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y, 
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);
void bicubic_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);
void cubicherm_2d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double ***f, double *value1, double *value2, double ***dfdx, double ***dfdy, double ***d2fdxdy);

// dfdx dfdy d2fdxdy can be supported or just NULL pointer.
// here 1d means for ONE 2d-functions. E.g f(x,y)
void bilenar_1d(double target_x, double target_y, int nx, double *x,  int ny, double *y, 
                double **f, double *value, double **dfdx, double **dfdy, double **d2fdxdy);

void cubicherm_1d(double target_x, double target_y, int nx, double *x,  int ny, double *y,
                double **f, double *value, double **dfdx, double **dfdy, double **d2fdxdy);


// void cubicherm_D1(double target_x, double target_y, int nx, double *x,  int ny, double *y,
//                 double **f, double *value, double **dfdx, double **dfdy, double **d2fdxdy);


// for 1D interpolation, x and f(x), we change to an abstrac interface.
// This is for 1D situation.
typedef void (*interp1d_eval_fun)(const void* interp_data, double x, double* value);
typedef void (*interp1d_deriv_fun)(const void* interp_data, double x, double* value);
typedef void (*interp1d_free_data)(void** ptr);


typedef struct 
{
  interp1d_eval_fun eval;
  interp1d_deriv_fun deriv;
  void* data;
  interp1d_free_data free_data;
  const char* name;
}Interp1DFunction;

typedef struct
{
  double* x;
  double* fx;
  double* dfdx;
  int nx;
}CubicHerm1dData;

void free_interp1d_function(Interp1DFunction* func);

//for cubicherm1D
void cubicherm1D_eval(const void* interp_data, double x, double* value);
void cubicherm1D_deriv(const void* interp_data, double x, double* value); 
Interp1DFunction* create_cubicherm1D_interp(double* x, double* fx, double* dfdx, int nx);
//void free_cubicherm1D_interp(Interp1DFunction* interp);
void cubicherm1D_free_data(void** ptr);
void modify_cubicherm1D_data(CubicHerm1dData* data, double* x, double* fx, double* dfdx, int nx);


/**************************************************
* To do, abstract for bilenar_1d, cubicherm_1d
***************************************************/
// for 2D interpolation, x,y and f(x,y), we change to an abstrac interface.
// This is for 2D situation.
typedef void (*interp2d_eval_fun)(const void* interp_data, double x, double y,double* value);
typedef void (*interp2d_deriv_fun)(const void* interp_data, double x, double y, double* value);
typedef void (*interp2d_free_data)(void** ptr);
typedef struct 
{
  interp2d_eval_fun eval;
  interp1d_deriv_fun deriv;
  void* data;
  interp1d_free_data free_data;
  const char* name;
}Interp2DFunction;

/**************************************************
* To do, abstract for bilenar_2d, cubicherm_2d
***************************************************/
// for 2D verctor interpolation vector, x,y and u(x,y) v(x,y), we change to an abstrac interface.
// This is for 2D situation.
typedef void (*interp2d_vec_eval_fun)(const void* interp_data, double x, double y, double* u, double* v);
typedef void (*interp2d_vec_deriv_fun)(const void* interp_data, double x, double y, double* u, double* v);
typedef void (*interp2d_vec_free_data)(void** ptr);
typedef struct 
{
  interp2d_vec_eval_fun eval;
  interp2d_vec_deriv_fun deriv;
  void* data;
  interp2d_vec_free_data free_data;
  const char* name;
}Interp2DFunction3D;

void Newton_Raphson_Method(double* x, const double* fx_target, Interp1DFunction* interp);
void Secant_Method(double* x, const double* fx_target, Interp1DFunction* interp);


#endif