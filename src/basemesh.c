#include"basemesh.h"

#ifndef EPS
#define EPS 1.0E-6
#endif

#ifndef NRELAX
#define NRELAX 500
#endif

#ifndef RLCEPT
#define RLCEPT 1.0E-6
#endif

#ifndef RELAX
#define RELAX 0.1
#endif


static double calc_length(double p1[2], double p2[2]) 
{
    return sqrt(pow(p1[0] - p2[0], 2.0) + pow(p1[1] - p2[1], 2.0));
}

static double dot_product(double a[2], double b[2]) 
{
    return a[0] * b[0] + a[1] * b[1];
}

static double cosine_term(double base[2], double point1[2], double point2[2]) 
{
    double vector1[2] = {point1[0] - base[0], point1[1] - base[1]};
    double vector2[2] = {point2[0] - base[0], point2[1] - base[1]};
    return dot_product(vector1, vector2);
}

void calc_ortho_CARRE(size_t n_point, double *length_prev_points, double **prev_point_coord,
                      double guard_top, double guard_end, double pasmin,
                      double *length_points, double **point_coord,
                      CarreOrthoValue *ortho_value)
{
  const double zero = 0.0;
  const double un = 1.0;

/********************************************
 * 
 * i-1 \        i+1/ previous curve
 *      \         /
 *       ----i----
 *           |
 *    -------I--------
 *   /                \
 *  / I-1           I+1\   cureve
 *
 *********************************************/
  double l1p, l2p, l12t, l12;
  int i = 0;
  l1p = calc_length(prev_point_coord[0], prev_point_coord[1]);

  l2p = calc_length(point_coord[0], point_coord[1]);

  l12t = calc_length(prev_point_coord[0], point_coord[0]);


  double g1 = guard_top;
  double g2 = guard_end;

  for (int i=1; i<n_point-1; i++)
  {
    double fac1, fac2, fac;
    double l12;
    if(guard_top > 0.0 && guard_end > 0.0)
    {
        fac1 = pow((g1/(g1 + length_prev_points[i])),2.0);
        fac2 = pow((g2/(g2 + length_prev_points[n_point-1] 
                             - length_prev_points[i])),2.0);
        fac = fac1*(un-fac2) + fac2;
        l12 = fac * length_prev_points[n_point-1] + (1.0-fac)*l12t;
    }
    else
    {
        l12=l12t;
    }

    double l1m = l1p;
    double l2m = l2p;

    l1p = calc_length(prev_point_coord[i], prev_point_coord[i+1]);
    l2p = calc_length(point_coord[i], point_coord[i+1]);
        
    double cs1, cs2, cs3, cs4;
    if (i - 1 >= 0 && i + 1 < n_point) 
    {
      cs1 = cosine_term(point_coord[i],prev_point_coord[i], point_coord[i-1]);
      cs2 = cosine_term(point_coord[i], prev_point_coord[i],point_coord[i+1]);
      cs3 = cosine_term(prev_point_coord[i], point_coord[i],point_coord[i-1]);
      cs4 = cosine_term(prev_point_coord[i], point_coord[i],point_coord[i+1]);
    }
    else
    {
      printf("Warning: in calc_ortho_CARRE, Out of bounds access at i=%d\n", i);
    }
        
    cs1=cs1/(l2m*l12);
    cs2=cs2/(l2p*l12);
    cs3=cs3/(l1m*l12);
    cs4=cs4/(l1p*l12);

    double f1,f2,f3;
    f1 = cs2+cs3-cs1-cs4;
    //this is according to the CARRE paper.
    // f2 = -(pow(g1/length_points[i],2.0) 
    //        +pow(g2/(length_points[n_point-1] - length_points[i]),2.0))
    //      *(length_points[i] 
    //        -length_prev_points[i]/length_prev_points[n_point-1]*length_points[n_point-1])
    //       /(pasmin + g1 + g2);
    //this is according to the CARRE source code clort.F
    f2 = -(pow(g1/length_prev_points[i],2.0) 
           +pow(g2/(length_prev_points[n_point-1] - length_prev_points[i]),2.0))
         *(length_points[i] 
           -length_prev_points[i]/length_prev_points[n_point-1]*length_points[n_point-1])
          /(pasmin + g1 + g2);
    f3 = (pow(pasmin/(length_points[i]-length_points[i-1]),2.0)
          -pow(pasmin/(length_points[i+1]-length_points[i]),2.0));
    
    ortho_value->ort[i] = f1 + f2 + f3;
    ortho_value->ortpur[i] = f1;
    ortho_value->propo[i] = f2;
    ortho_value->varr[i] = f3;
    ortho_value->tot[i] = f1 + f2 + f3;
  }
}


void calc_points_CARRE(CarreMeshTube *tube)
{
  CarreOrthoValue tmp_ortho_value, ortho_value;
  allocate_CarreOrthoValue(tube->n_point, &tmp_ortho_value);
  allocate_CarreOrthoValue(tube->n_point, &ortho_value);
  size_t n_point = tube->n_point;

  if(NRELAX >= 0)
  {
    double d1 = 0.0;
    size_t ipol1 = 1;
    size_t ipoln = n_point - 1;
//1. we first arrange the points proportionally to those of the previous line.
    //tmp_length_points is coresponding to l1 in mailrg.F
    //tube->length_points is coresponding to l2 mailrg.F
    double *tmp_length_points = (double *)calloc(n_point, sizeof(double));
    double **tmp_point_coord = allocate_2d_array(n_point,2);

    tmp_length_points[0] = 0.0; 
    double length = long_CARRE(tube->curve, tube->n_curve);
    tmp_length_points[n_point-1] = length;
    double prev_length = tube->length_prev_points[n_point-1];
    
    // the first and the last point of mesh points are the same with curve
    tube->point_coord[0][0] = tube->curve[0][0];
    tube->point_coord[0][1] = tube->curve[0][1];
    tube->point_coord[n_point-1][0] = tube->curve[tube->n_curve-1][0];
    tube->point_coord[n_point-1][1] = tube->curve[tube->n_curve-1][1];

    tmp_point_coord[0][0] = tube->curve[0][0];
    tmp_point_coord[0][1] = tube->curve[0][1];
    tmp_point_coord[n_point-1][0] = tube->curve[tube->n_curve-1][0];
    tmp_point_coord[n_point-1][1] = tube->curve[tube->n_curve-1][1];
    printf("debug npoint-1: %d, n_curve-1: %d \n", (int)(n_point-1), (int)(tube->n_curve-1));

//    printf("debug in calc_points_CARRE line 149\n");

    // the remaing point from 1 to n_point-2
    for(int ipol = ipol1; ipol<ipoln; ipol++)
    {
      d1 = ruban_CARRE(tube->prev_curve, tube->n_prev_curve, 
                       tube->prev_point_coord[ipol],d1);
      tmp_length_points[ipol] = (tube->length_prev_points[ipol]/prev_length)*length;
 
      coord_CARRE(tube->curve, tube->n_curve, tmp_length_points[ipol], tmp_point_coord[ipol]);
    }
    //debug******************************************
    const char *ini_point_name="ini_point";
    FILE* ini_point = fopen(ini_point_name, "w");
    for (int i=0; i<n_point; i++)
    {
      fprintf(ini_point, "%.12f %.12f\n", tmp_point_coord[i][0],tmp_point_coord[i][1]);
    }
    fclose(ini_point);
    printf("write the tracing line in %s\n", ini_point_name);
    //***********************************************************************

 //   printf("debug in calc_points_CARRE line 159\n");

    if(NRELAX > 0)
    {
//2. we initialize the function which must be zero for an orthogonal distribution.
      calc_ortho_CARRE(tube->n_point,tube->length_prev_points, tube->prev_point_coord,
                      tube->guard_top, tube->guard_end, tube->pasmin,
                      tmp_length_points,tmp_point_coord,
                      &tmp_ortho_value);
      printf("debug origin ort: %.10f\n",tmp_ortho_value.ort[n_point-1]);

//3. we proceed to a first displacement of the nodes

      tube->length_points[0] = 0;
      tube->length_points[n_point-1] = length;
  
      for(int ipol = 1; ipol<ipoln; ipol++)
      {
        if(tmp_ortho_value.ort[ipol] > 0.0)
        {
          tube->length_points[ipol] = 0.9*tmp_length_points[ipol] + 0.1*tmp_length_points[ipol+1];
        }
        else
        {
          tube->length_points[ipol] = 0.9*tmp_length_points[ipol] + 0.1*tmp_length_points[ipol-1];
        }
        coord_CARRE(tube->curve,tube->n_curve,tube->length_points[ipol], tube->point_coord[ipol]);
      }

    //debug******************************************
    const char *firstname="first_replace";
    FILE* first = fopen(firstname, "w");
    for (int i=0; i<n_point; i++)
    {
      fprintf(first, "%.12f %.12f\n", tube->point_coord[i][0], tube->point_coord[i][1]);
    }
    fclose(first);
    printf("write the tracing line in %s\n", firstname);
    //************************************************************
    //todo: store the ortho value for the whole region
    // somort(ir)= somort(ir)+ (ort1(ipol)/nppol)
    // somortpur(ir)= somortpur (ir)+ (ortpur1(ipol)/nppol)
    // sompropo(ir)= sompropo(ir)+ (propo1(ipol)/nppol)
    // somvarr(ir)= somvarr(ir)+ (varr1(ipol)/nppol)
    // somtot(ir)= somtot(ir)+ (tot1(ipol)/nppol)
    
//4. we relax the points iteratively to achieve the best possible orthogonality
      double ortmax = 0.0;
      
      for (int i = 0; i<NRELAX; i++)
      {
        calc_ortho_CARRE(tube->n_point,tube->length_prev_points, tube->prev_point_coord,
                         tube->guard_top, tube->guard_end, tube->pasmin,
                         tube->length_points, tube->point_coord,
                         &ortho_value);
        ortmax = 0.0;

        for (int ipol=ipol1; ipol<ipoln; ipol++)
        {
          double ortho_current = ortho_value.ort[ipol];
          double ortho_diff = ortho_current - tmp_ortho_value.ort[ipol];
          double length_diff = tube->length_points[ipol] - tmp_length_points[ipol];
          double del = 0.0;
          double ortho_diff_abs = fabs(ortho_diff);

          if (fabs(ortho_current) <= RLCEPT)
          {
            continue;
          }

          if(ortho_diff_abs > RLCEPT*RLCEPT)
          {
            del = -ortho_current*length_diff/ortho_diff;
          }
          else
          {
            del = 0.0;
          } 

          if(del > 0.0) 
          {
            del = min(del, RELAX*(tube->length_points[ipol+1]-tube->length_points[ipol]));
          }
          else
          {
            del = max(del, RELAX*(tube->length_points[ipol-1]-tube->length_points[ipol]));
          }

          if(fabs(del) > RLCEPT*RLCEPT)
          {
            tmp_length_points[ipol] = tube->length_points[ipol];
            tmp_ortho_value.ort[ipol] = ortho_value.ort[ipol];
            tube->length_points[ipol] = tmp_length_points[ipol]+del;
          }
          // printf("debug in calc_points_CARRE line 241\n"); 
          //printf("del: %.10f\n", del);
          //printf("ipol: %d, length_point %.10f\n",ipol, tube->length_points[ipol]);
          coord_CARRE(tube->curve,tube->n_curve,tube->length_points[ipol], tube->point_coord[ipol]);
          ortmax=max(ortmax,fabs(ortho_current));
        }
        printf("debug: ortmax %.10f\n",ortmax);
        printf("i: %d\n",i);
        if ((ortmax < RLCEPT))
        {
          printf("Finish optimized\n");
          break;
        }
      }
      if (ortmax > RLCEPT)
      {
        printf("WARING: After optimization, the orthognoality is not good\n");
        printf("please adjust the parameters\n");
      }
      //TODO write the orthognonal values
      //   do ipol=ipol1,ipoln
      //   somortp(ir)= somortp(ir)+ (ort2(ipol)/nppol)
      //   somortpurp(ir)= somortpurp(ir)+ (ortpur2(ipol)/nppol)
      //   sompropop(ir)= sompropop(ir)+ (propo2(ipol)/nppol)
      //   somvarrp(ir)= somvarrp(ir)+ (varr2(ipol)/nppol)
      //   somtotp(ir)= somtotp(ir)+ (tot2(ipol)/nppol)
      //   enddo
      free(tmp_length_points);
      free_2d_array(tmp_point_coord);
    }
    else
    {
      printf("not support NRELAX=0\n");
    }
  }
  free_CarreOrthoValue(&tmp_ortho_value);
  free_CarreOrthoValue(&ortho_value);

}

int indcrb_CARRE(double **curve, size_t n_curve, double point[2], double d)
{
  double mumin = 1.0;
  int indcrb = 0;
  double ax, ay, bx, by, dist, mu;
  double dd = 0.0; 

  ax = curve[0][0] - point[0];
  ay = curve[0][1] - point[1];
  dist = sqrt(ax * ax + ay * ay);

  double first_segment_length = sqrt(pow(curve[1][0] - curve[0][0], 2) 
                                     +pow(curve[1][1] - curve[0][1], 2));
  if (dist < EPS && first_segment_length > d) 
  {
    return 1;
  }

  for (size_t i = 0; i < n_curve - 1; i++) 
  {
    dd += sqrt(pow(curve[i + 1][0] - curve[i][0], 2) 
               +pow(curve[i + 1][1] - curve[i][1], 2));

    if (dd >= d)
    {
      ax = curve[i][0] - point[0];
      ay = curve[i][1] - point[1];
      bx = curve[i + 1][0] - point[0];
      by = curve[i + 1][1] - point[1];

      dist = sqrt(bx * bx + by * by);
      if (dist < EPS) 
      {
        return (i < n_curve - 2) ? i+1 : i; 
      }
      else
      {
        double norm_a = sqrt(ax * ax + ay * ay);
        double norm_b = sqrt(bx * bx + by * by);
        mu = (ax * bx + ay * by) / (norm_a * norm_b);     
        if (mu < mumin) 
        {
          mumin = mu;
          indcrb = i;
        }
      }
    } 
  }
  return indcrb;    
}

double long_CARRE(double **curve, size_t n) 
{
    double total_length = 0.0; 
    double dist;
    for (int i = 0; i < n - 1; i++) 
    {
    //  printf("debug i: %d\n", i);  
      dist = sqrt(pow(curve[i + 1][0] - curve[i][0], 2.0) + pow(curve[i + 1][1] - curve[i][1], 2.0));
      total_length += dist;
    }

    return total_length;
}

double ruban_CARRE(double **curve, size_t n, double point[2], double d) 
{
    double ruban = 0.0; 
    double dist = 0.0;  
    int ind;           
    double dx, dy;      

    ind = indcrb_CARRE(curve, n, point, d);

    dist = long_CARRE(curve, ind);

    dx = point[0] - curve[ind][0];
    dy = point[1] - curve[ind][1];

    ruban = dist + sqrt(dx * dx + dy * dy);

    return ruban;
}

void coord_CARRE(double **curve, size_t n, double d, double point[2])
{
  double dist = 0.0;
  double distot = 0.0;
  point[0] = 0.0;
  point[1] = 0.0;

  for (size_t i = 0; i < n - 1; i++) 
  {
    double dx = curve[i + 1][0] - curve[i][0];
    double dy = curve[i + 1][1] - curve[i][1];
    dist = sqrt(dx * dx + dy * dy);
    distot += dist;
    if (distot >= d) 
      {
        double ratio = (distot - d) / dist;
        point[0] = curve[i + 1][0] - ratio * dx;
        point[1] = curve[i + 1][1] - ratio * dy;
        return;
      }
  }

  if (fabs(distot - d) < EPS)
  {
    point[0] = curve[n - 1][0];
    point[1] = curve[n - 1][1];
    return;
  }

  printf("Error: Required distance (%.6f) exceeds total curve length (%.6f)\n", d, distot);
  exit(1);

}

void allocate_CarreOrthoValue(int n, CarreOrthoValue *ortho_value)
{
    ortho_value->n = n;
    ortho_value->ort = (double *)calloc(n, sizeof(double));
    ortho_value->ortpur = (double *)calloc(n, sizeof(double));
    ortho_value->propo  = (double *)calloc(n, sizeof(double));
    ortho_value->varr = (double *)calloc(n, sizeof(double));
    ortho_value->tot = (double *)calloc(n, sizeof(double));
}

void free_CarreOrthoValue(CarreOrthoValue *ortho_value)
{
    free(ortho_value->ort) ;
    free(ortho_value->ortpur); 
    free(ortho_value->propo); 
    free(ortho_value->varr);
    free(ortho_value->tot);
}