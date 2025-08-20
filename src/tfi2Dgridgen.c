#include "tfi2Dgridgen.h"

#define DISMIN 0.05
#define MAX_LAPLACE_ITER 5000

/* check whether the point (x1,y1) is identical with (x2,y2)
 * return 1 if they are equal, 0 otherwise
 */
static int points_equal(double x1, double y1, double x2, double y2)
{
  if (fabs(x1 - x2) < 1.0e-15 && fabs(y1 - y2) < 1.0e-15) 
  {
    return 1;
  } 
  else 
  {
    return 0;
  }
}

void check_tfi_boundary(Curve* cur_l, Curve* cur_r, Curve* cur_b, Curve* cur_t)
{
  size_t nL = get_curve_npt(cur_l);
  size_t nR = get_curve_npt(cur_r);
  size_t nB = get_curve_npt(cur_b);
  size_t nT = get_curve_npt(cur_t);

  if (nL == 0 || nR == 0 || nB == 0 || nT == 0) 
  {
    fprintf(stderr, "One or more boundary curves have zero points.\n");
    exit(EXIT_FAILURE);
  }

  /* Left first vs Bottom first */
  if (!points_equal(get_curve_x(cur_l, 0), get_curve_y(cur_l, 0),
               get_curve_x(cur_b, 0), get_curve_y(cur_b, 0))) 
  {
    fprintf(stderr, "The first point of Left boundary is not identical with the first of Bottom boundary.\n");
    exit(EXIT_FAILURE);
  }

  /* Bottom last vs Right first */
  if (!points_equal(get_curve_x(cur_b, nB - 1), get_curve_y(cur_b, nB - 1),
                 get_curve_x(cur_r, 0),      get_curve_y(cur_r, 0))) 
  {
    fprintf(stderr, "The last point of Bottom boundary is not identical with the first of Right boundary.\n");
    exit(EXIT_FAILURE);
  }

  /* Left last vs Top first */
  if (!points_equal(get_curve_x(cur_l, nL - 1), get_curve_y(cur_l, nL - 1),
               get_curve_x(cur_t, 0),      get_curve_y(cur_t, 0))) 
  {
    fprintf(stderr, "The last point of Left boundary is not identical with the first of Top boundary.\n");
    exit(EXIT_FAILURE);
  }

  /* Top last vs Right last */
  if (!points_equal(get_curve_x(cur_t, nT - 1), get_curve_y(cur_t, nT - 1),
               get_curve_x(cur_r, nR - 1), get_curve_y(cur_r, nR - 1))) 
  {
    fprintf(stderr, "The last point of Top boundary is not identical with the last of Right boundary.\n");
    exit(EXIT_FAILURE);
  }
}


void generate_2Dgrid_default_TFI(TwoDimGrid* grid,
                                 Curve* cur_b, Curve* cur_t, Curve* cur_l, Curve* cur_r,
                                 const double* distb_cur_b, int nbottom,
                                 const double* distb_cur_t, int ntop,
                                 const double* distb_cur_l, int nleft, 
                                 const double* distb_cur_r, int nright)
{
  if(!grid || !cur_l || !cur_r || !cur_b || !cur_t || !distb_cur_l || !distb_cur_r || !distb_cur_b || !distb_cur_t )
  {
    fprintf(stderr, "Empty inputs for generate_2Dgrid_default_TFI.\n");
    exit(EXIT_FAILURE);
  }

  if(nleft!=nright || nbottom!=ntop)
  {
    fprintf(stderr, "The grid points are not equal at the boundaries.\n");
    exit(EXIT_FAILURE);
  }

  if(nbottom!=grid->npol || nleft!=grid->nrad)
  {
    fprintf(stderr, "Thes siez of 2D grid are not consistant with inputs in generate_2Dgrid_default_TFI.\n");
    exit(EXIT_FAILURE);
  }

  if(fabs(distb_cur_l[0])>EPSILON_15 ||
     fabs(distb_cur_r[0])>EPSILON_15 ||
     fabs(distb_cur_b[0])>EPSILON_15 ||
     fabs(distb_cur_t[0])>EPSILON_15
     )
  {
    fprintf(stderr, "The normalized distrions are not from ZERO.\n");
    exit(EXIT_FAILURE);
  }

  if(fabs(distb_cur_l[nleft-1]-1.0)>EPSILON_15 ||
     fabs(distb_cur_r[nright-1]-1.0)>EPSILON_15 ||
     fabs(distb_cur_b[nbottom-1]-1.0)>EPSILON_15 ||
     fabs(distb_cur_t[ntop-1]-1.0)>EPSILON_15
     )
  {
    fprintf(stderr, "The normalized distrions are not end with ONE.\n");
    exit(EXIT_FAILURE);
  }

  check_tfi_boundary(cur_l, cur_r, cur_b, cur_t);
  
  #ifdef DEBUG
    printf("The inputs for generate_2Dgrid_default_TFI are ready.\n");
  #endif
// Begin to construct the 2d grid which includes the normalized distribution u and v

  TwoDimGrid* norm_grid=create_2Dgrid_poloidal_major(nbottom, nleft);

  for (size_t i = 0; i < nbottom; i++)
  {
    set_point_2Dgrid(norm_grid, i,       0, distb_cur_b[i], 0.0);
    set_point_2Dgrid(norm_grid, i, nleft-1, distb_cur_t[i], 1.0);
  }
  
  for (size_t j = 0; j < nleft; j++)
  {
    set_point_2Dgrid(norm_grid,         0, j, 0.0, distb_cur_l[j]);
    set_point_2Dgrid(norm_grid, nbottom-1, j, 1.0, distb_cur_r[j]);
  }
 
  
// BoundaryBlendedControlFunction refer to https://github.com/amalinadhi/2DTFIGridGeneration/blob/master/core/CoreMethod.py
  int iMax=nbottom;
  int jMax=nleft;
  
  for(int j=1; j<jMax-1; j++)  
  {
    for(int i=1; i<iMax-1; i++)
    {
      double part1 = (1.0 - get_y_2Dgrid(norm_grid,0,j)) * get_x_2Dgrid(norm_grid,i,0)
                   + get_y_2Dgrid(norm_grid,0,j) * get_x_2Dgrid(norm_grid,i,jMax-1);
        
      double part2 = 1.0 - (get_x_2Dgrid(norm_grid,i,jMax-1) - get_x_2Dgrid(norm_grid,i,0)) *
                           (get_y_2Dgrid(norm_grid,iMax-1,j) - get_y_2Dgrid(norm_grid,0,j));
        
      double part3 = (1.0 - get_x_2Dgrid(norm_grid,i,0)) * get_y_2Dgrid(norm_grid,0,j)
                   + get_x_2Dgrid(norm_grid,i,0) * get_y_2Dgrid(norm_grid,iMax-1,j);
        
      double part4 = 1.0 - (get_y_2Dgrid(norm_grid,iMax-1,j) - get_y_2Dgrid(norm_grid,0,j)) *
                           (get_x_2Dgrid(norm_grid,i,jMax-1) - get_x_2Dgrid(norm_grid,i,0));
        
      set_point_2Dgrid(norm_grid, i, j, part1/part2, part3/part4);
    }
  }
  #ifdef DEBUG
    write_2Dgrid(norm_grid, "NORM_DISTRIBUTION_TFI");
  #endif
  // Calculate the grid point in the four boudanries according to the normalized distributions and curves
  double tot_len = total_length_curve(cur_b);
  for(int i=0;i<iMax;i++)
  {
    double len = tot_len*distb_cur_b[i];
    CurvePoint point;
    coordnates_in_curve(cur_b, len, &point);
    set_point_2Dgrid(grid, i, 0, point.x, point.y);
  }

  tot_len=total_length_curve(cur_t);
  for(int i=0;i<iMax;i++)
  {
    double len = tot_len*distb_cur_t[i];
    CurvePoint point;
    coordnates_in_curve(cur_t, len, &point);
    set_point_2Dgrid(grid, i, jMax-1, point.x, point.y);
  }

  tot_len = total_length_curve(cur_l);
  for(int j=0;j<jMax;j++)
  {
    double len = tot_len*distb_cur_l[j];
    CurvePoint point;
    coordnates_in_curve(cur_l, len, &point);
    set_point_2Dgrid(grid, 0, j, point.x, point.y);
  }

  tot_len = total_length_curve(cur_r);
  for(int j=0;j<jMax;j++)
  {
    double len = tot_len*distb_cur_r[j];
    CurvePoint point;
    coordnates_in_curve(cur_r, len, &point);
    set_point_2Dgrid(grid, iMax-1, j, point.x, point.y);
  }
 
  // Transfinite Interpolation for the 2d grid points from [1:iMax-1] and [1:jMax-1] 0-based
  for(int j=1;j<jMax-1;j++)
  {
    for(int i=1;i<iMax-1;i++)
    {
      double u = get_x_2Dgrid(norm_grid, i, j);
      double v = get_y_2Dgrid(norm_grid, i, j);

      double Ux = (1-u) * get_x_2Dgrid(grid, 0, j) 
                +    u  * get_x_2Dgrid(grid, iMax-1, j);
        
      double Vx = (1-v) * get_x_2Dgrid(grid, i, 0)
                +    v  * get_x_2Dgrid(grid, i, jMax-1);
        
      double UVx = u*v         * get_x_2Dgrid(grid, iMax-1, jMax-1)
                 + u*(1-v)     * get_x_2Dgrid(grid, iMax-1, 0)
                 + (1-u)*v     * get_x_2Dgrid(grid, 0, jMax-1)
                 + (1-u)*(1-v) * get_x_2Dgrid(grid, 0, 0);

      double Uy = (1-u) * get_y_2Dgrid(grid, 0, j)
                +    u  * get_y_2Dgrid(grid, iMax-1, j);
        
      double Vy = (1-v) * get_y_2Dgrid(grid, i, 0)
                +    v  * get_y_2Dgrid(grid, i, jMax-1);
        
      double UVy = u*v         * get_y_2Dgrid(grid, iMax-1, jMax-1)
                 + u*(1-v)     * get_y_2Dgrid(grid, iMax-1, 0)
                 + (1-u)*v     * get_y_2Dgrid(grid, 0, jMax-1)
                 + (1-u)*(1-v) * get_y_2Dgrid(grid, 0, 0);

      set_point_2Dgrid(grid,i,j,Ux + Vx - UVx, Uy + Vy - UVy);
    }
  }

  #ifdef DEBUG
    printf("Finish the grid generation by TFI.\n");
  #endif

  free_2Dgrid(norm_grid);
}

// void fix_out_of_bounds_2Dgrid(TwoDimGrid* grid)
// {


// }

static int is_point_on_edge(double x, double y, CurvePoint p1, CurvePoint p2) 
{

    double min_x = fmin(p1.x, p2.x);
    double max_x = fmax(p1.x, p2.x);
    double min_y = fmin(p1.y, p2.y);
    double max_y = fmax(p1.y, p2.y);
    
    if (x < min_x - EPSILON_12 || x > max_x + EPSILON_12 || 
        y < min_y - EPSILON_12 || y > max_y + EPSILON_12) 
    {
      return 0;
    }
    
    double cross_product = (y - p1.y) * (p2.x - p1.x) - (x - p1.x) * (p2.y - p1.y);
    return fabs(cross_product) <= EPSILON_12;
}

//check whether a point is inside a closed region which have ncurve curves or on the edge
//if true, return 1, else retrun 0
static int is_point_inside(double x, double y, Curve** boundary, int n_boundary)
{
  int crossings = 0;
  for(int c=0;c<n_boundary;c++)
  {
    int n_point  =boundary[c]->n_point;
    for (int i=0;i<n_point-1;i++)
    {
      CurvePoint p1=boundary[c]->points[i];
      CurvePoint p2=boundary[c]->points[i+1];
      
      //check whether is on the edge
      if (is_point_on_edge(x, y, p1, p2)) 
      {
        return 1; 
      }

      if (p1.y > p2.y) 
      {
        CurvePoint tmp = p1;
        p1 = p2;
        p2 = tmp;
      }
      if (y > p1.y && (y < p2.y || fabs(y - p2.y) <= EPSILON_12)) 
      {
        if (fabs(p2.y - p1.y) > EPSILON_12) 
        {
          double x_inter = p1.x + (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
          if (x < x_inter) 
          {
            crossings++;
          }
        }
      }
    }
  }
  return crossings%2;
}

int count_points_inside(TwoDimGrid* grid, Curve** boundary, int n_boundary)
{
  int np=grid->npol;
  int nr=grid->nrad;
  int number = 0;
  for(int j=1;j<nr-1;j++)
  {
    for(int i=1;i<np-1;i++)
    {
      if(is_point_inside(get_x_2Dgrid(grid,i,j),get_y_2Dgrid(grid,i,j),boundary,n_boundary)==0)
      {
        number++;
      }
    }
  }
  return number;
}

int try_neighbor_averaging(TwoDimGrid* grid, Curve** boundary, int n_boundary, int max_iter,double drmin)
{
  for (int iter = 0; iter < max_iter; iter++) 
  {
    int np=grid->npol;
    int nr=grid->nrad;
      
    for (int i = 1; i < np-1; i++) 
    {
      for (int j = 1; j < nr-1; j++) 
      {
        double x=get_x_2Dgrid(grid,i,j)-get_x_2Dgrid(grid,i,j-1);
        double y=get_y_2Dgrid(grid,i,j)-get_y_2Dgrid(grid,i,j-1);
        if (!is_point_inside(get_x_2Dgrid(grid,i,j), get_y_2Dgrid(grid,i,j), boundary, n_boundary)
            || hypot(x,y)<drmin) 
        {
          for(int k=j;k<nr-1;k++)
          {
            double x_avg = 0.25 * (get_x_2Dgrid(grid,i+1,k) + get_x_2Dgrid(grid,i-1,k) + 
                                   get_x_2Dgrid(grid,i,k+1) + get_x_2Dgrid(grid,i,k-1));
            double y_avg = 0.25 * (get_y_2Dgrid(grid,i+1,k) + get_y_2Dgrid(grid,i-1,k) +
                                   get_y_2Dgrid(grid,i,k+1) + get_y_2Dgrid(grid,i,k-1));
            set_point_2Dgrid(grid,i,k,x_avg,y_avg);
          }
        break;
        }
      }
    }
  }
  return count_points_inside(grid, boundary, n_boundary);
}




void optimized_neu_2Dgrid(TwoDimGrid* grid)
{
  Curve* boundary[4];
  for(int i=0;i<4;i++)
  {
    boundary[i]=create_curve(0);
  }
  int np=grid->npol;
  int nr=grid->nrad;

  for(int i=0;i<np;i++)
  {
    add_last_point_curve(boundary[0],get_x_2Dgrid(grid,i,0),get_y_2Dgrid(grid,i,0));
    add_last_point_curve(boundary[1],get_x_2Dgrid(grid,i,nr-1),get_y_2Dgrid(grid,i,nr-1));
  }

  for(int j=0;j<nr;j++)
  {
    add_last_point_curve(boundary[2],get_x_2Dgrid(grid,0,j),get_y_2Dgrid(grid,0,j));
    add_last_point_curve(boundary[3],get_x_2Dgrid(grid,np-1,j),get_y_2Dgrid(grid,np-1,j));
  }

  // for(int i=0;i<10;i++)
  // {
  //   try_avoid_drmin(grid, 0.005);
  // }
  #ifdef DEBUG
    write_curve("neu_2D_boudnar_0",boundary[0]);
    write_curve("neu_2D_boudnar_1",boundary[1]);
    write_curve("neu_2D_boudnar_2",boundary[2]);
    write_curve("neu_2D_boudnar_3",boundary[3]);
  #endif
  if(count_points_inside(grid,boundary,4)==0)
  {
    printf("All the points are inside the boundaries and no need to optimized.\n");
  }
  else if (try_neighbor_averaging(grid,boundary, 4, 50, DISMIN)==0)
  {
    printf("After neighbor_averaging, all the points are inside the boundaries.\n");
    LaplaceSmoothing(grid, 1.2, 1.0E-3);
  }
  else
  {
    fprintf(stderr, "Still have points are outside the boundaries.\n");
    fprintf(stderr, "Please adjust top boundary! In the future, an automatic optimization will come.\n");
    exit(EXIT_FAILURE);
  }

  for(int i=0;i<4;i++)
  {
    free_curve(boundary[i]);
  }
}

int LaplaceSmoothing(TwoDimGrid* grid, double omega, double targetError)
{
  int iteration=0;
  double error = 1.0;
  double lastRValue =1.0;
  double residual;

  int np = grid->npol;
  int nr = grid->nrad;
  while(error>targetError)
  {
    double Rx_sum = 0.0;
    double Ry_sum = 0.0;
    for(int j=1;j<nr-1;j++)
    {
      for(int i=1;i<np-1;i++)
      {
        // Calculate grid metric parameters
        double xXi = (get_x_2Dgrid(grid,i+1,j) - get_x_2Dgrid(grid,i-1,j))/2;
        double yXi = (get_y_2Dgrid(grid,i+1,j) - get_y_2Dgrid(grid,i-1,j))/2;
        double xEta = (get_x_2Dgrid(grid,i,j+1) - get_x_2Dgrid(grid,i,j-1))/2;
        double yEta = (get_y_2Dgrid(grid,i,j+1) - get_y_2Dgrid(grid,i,j-1))/2;

        // Calculate Jacobian and metric coefficients
        double J = xXi*yEta - xEta*yXi;
        double alpha = xEta*xEta + yEta*yEta;
        double beta = xXi*xEta + yXi*yEta;
        double gamma = xXi*xXi + yXi*yXi;

        // Calculate X-direction residual
        double Rx1 = alpha*(get_x_2Dgrid(grid,i+1,j) - 2*get_x_2Dgrid(grid,i,j) + get_x_2Dgrid(grid,i-1,j));
        double Rx2 = (-0.5)*beta*(get_x_2Dgrid(grid,i+1,j+1) - get_x_2Dgrid(grid,i-1,j+1) - 
                                  get_x_2Dgrid(grid,i+1,j-1) + get_x_2Dgrid(grid,i-1,j-1));
        double Rx3 = gamma*(get_x_2Dgrid(grid,i,j+1) - 2*get_x_2Dgrid(grid,i,j) + get_x_2Dgrid(grid,i,j-1));
        double Rx = Rx1 + Rx2 + Rx3;

        // Calculate Y-direction residual
        double Ry1 = alpha*(get_y_2Dgrid(grid,i+1,j) - 2*get_y_2Dgrid(grid,i,j) + get_y_2Dgrid(grid,i-1,j));
        double Ry2 = (-0.5)*beta*(get_y_2Dgrid(grid,i+1,j+1) - get_y_2Dgrid(grid,i-1,j+1) - 
                                  get_y_2Dgrid(grid,i+1,j-1) + get_y_2Dgrid(grid,i-1,j-1));
        double Ry3 = gamma*(get_y_2Dgrid(grid,i,j+1) - 2*get_y_2Dgrid(grid,i,j) + get_y_2Dgrid(grid,i,j-1));
        double Ry = Ry1 + Ry2 + Ry3;

        // Accumulate residuals for convergence check
        Rx_sum += Rx;
        Ry_sum += Ry;

        // Update X and Y coordinates
        double x_new = get_x_2Dgrid(grid,i,j) + omega*(Rx/(2*(alpha + gamma)));
        double y_new = get_y_2Dgrid(grid,i,j) + omega*(Ry/(2*(alpha + gamma)));
        set_point_2Dgrid(grid, i, j, x_new, y_new);
      }
    }
    // Find residual
    double currentRValue = sqrt(Rx_sum*Rx_sum + Ry_sum*Ry_sum);
    error = fabs(lastRValue - currentRValue);
    iteration++;

    if (iteration > MAX_LAPLACE_ITER) 
    {
      printf("Reach the maximum iterations in LaplaceSmoothing.\n");
      break;
    }
    lastRValue = currentRValue;
  }
  printf("Finish LaplaceSmoothing with omega %lf, targetError %lf.\n", omega, targetError);
  return iteration;
}