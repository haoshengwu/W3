#include "xpoint.h"

//use to determine whethere the found X-point is near the estimated Xpoint within the
//range r: [R-DELTA, R+DELTA] z: [Z-DELTA, Z+DELTA]
#ifndef DELTA
#define DELTA 0.1
#endif

#ifndef LEARNINGRATE
#define LEARNINGRATE 5E-6
#endif

#ifndef MAXITER
#define MAXITER 10000
#endif

#ifndef TOLERANCE
#define TOLERANCE 8.0E-4
#endif

// private functions for find Xpoint
static int check_xpt_rectangular(const Equilibrium *equilib, const XPointInfo xpt);
static int check_xpt_est_range(const Equilibrium *equilib, int xpoint_number, double **est_xpoint_pos, const XPointInfo xpt);
static int check_xpt_levels(const Equilibrium *equilib,int cx1,int cy1,int cx2,int cy2,int x0,int y0,int bMinMax);
static void calculate_xpt_level(Equilibrium *equilib, int xpoint_number, double **est_xpoint_pos,
                                interpl_1D_f interpl_1D_f, interpl_2D_f interpl_2D_f, _XPointInfo *xpt_array);


void find_xpoint(Equilibrium *equilib, int xpoint_number, double **est_xpoint_pos, 
                 interpl_1D_f interpl_1D_f, interpl_2D_f interpl_2D_f, _XPointInfo *xpt_array)
{
  // est_xpoint_pos is the estimated postions of xpoint. est_xpoint_pos[:][0] is R, est_xpoint_pos[:][1] is Z. 
  // xpoint_number is the expected X-point number. Currentlu, only suppot ONE Xpoint. 
  // In the future, multiple Xpoints will be considered.
  // this algorithm is from DivGeo FindXPointRects and FindXPointCenter xpoint.h and xpoint.c files.
  int i, j, s;
  int ixpt=-1; //store the current X-pints numbers 
  XPointInfo xp, xpC; // xpC is for every potentional xpoint. xp is the one store the X-point.

  xpC = malloc(sizeof(*xpC));
//  xp = malloc(sizeof(*xp));

  for (s = 1; s < 3; s++)
  {
    for (i = 1; i < equilib->nw - 1 - s; i++)
    {
      for (j = 1; j < equilib->nh - 1 - s; j++)
      {
        xpC->cx1 = i,
        xpC->cy1 = j;
        xpC->cx2 = i + s;
        xpC->cy2 = j + s;

        if (check_xpt_rectangular(equilib, xpC))
        {
          continue;
        };

        //printf("Find a potential equ rectangular contain a Xpoint:\n");
        //printf("crx1: %i, crx2: %i, cry1: %i, cry2: %i\n", xpC->cx1, xpC->cx2, xpC->cy1, xpC->cy2);
        //printf("x1: %.12f y1: %.12f\n", equilib->r[xpC->cx1], equilib->z[xpC->cy1]);
        //printf("minMax[0:3]: %.12f %.12f %.12f %.12f\n", xpC->minMax[0].lvl, xpC->minMax[1].lvl,
        //                                                 xpC->minMax[2].lvl, xpC->minMax[3].lvl);
        // store the value in xp
        if (!check_xpt_est_range(equilib, xpoint_number, est_xpoint_pos, xpC))
        {
          if (ixpt == -1 || 
            ((ixpt != -1) && abs(xpC->cx1 - xpt_array[ixpt].cx1) > s && abs(xpC->cy1 - xpt_array[ixpt].cy1) > s)) 
            {
              ixpt++;
              xpt_array[ixpt].cx1 = xpC->cx1;
              xpt_array[ixpt].cx2 = xpC->cx2;
              xpt_array[ixpt].cy1 = xpC->cy1;
              xpt_array[ixpt].cy2 = xpC->cy2;
              xpt_array[ixpt].lvlMin = xpC->lvlMin;
              xpt_array[ixpt].lvlMax = xpC->lvlMax;
              xpt_array[ixpt].minMax[0] = xpC->minMax[0];
              xpt_array[ixpt].minMax[1] = xpC->minMax[1];
              xpt_array[ixpt].minMax[2] = xpC->minMax[2];
              xpt_array[ixpt].minMax[3] = xpC->minMax[3];
              printf("ixpt: %d\n",ixpt);
              printf("crx1: %i, crx2: %i, cry1: %i, cry2: %i\n", xpC->cx1, xpC->cx2, xpC->cy1, xpC->cy2);
              printf("x1: %.12f y1: %.12f\n", equilib->r[xpC->cx1], equilib->z[xpC->cy1]);
              printf("minMax[0:3]: %.12f %.12f %.12f %.12f\n", xpC->minMax[0].lvl, xpC->minMax[1].lvl,
                                                        xpC->minMax[2].lvl, xpC->minMax[3].lvl);
            }
          if (ixpt >= xpoint_number)
          {
            printf("Expected Xpoint number: %d\n", xpoint_number);
            printf("More xpionts are found. Please check the equilibriun/estimated xpoints positions\n");
            exit(EXIT_FAILURE);;
          }
        }
      }
    }
  }
  if (ixpt != xpoint_number-1)
  {
    printf("Expected Xpoint number: %d\n", xpoint_number);
    printf("Current Xpoint number: %d\n", ixpt+1);
    printf("Less xpionts are found. Please check the equilibriun/estimated xpoints positions\n");
    exit(1);
  }
  calculate_xpt_level(equilib, xpoint_number, est_xpoint_pos, interpl_1D_f, interpl_2D_f, xpt_array);
  printf("Finish find the X-points.\n");
  return;
}


// refer to DivGeo source code [static int CheckXPointRect(Equil eq,XPointTest xpt)]
// return 0 means find the xpoint in the rectangular, -1 means fail.
static int check_xpt_rectangular(const Equilibrium *equilib, const XPointInfo xpt)
{
  struct XPointExtremum p[6];
  int x, y, ox, oy, n, d, i;
  double lvl, lvl1; // the psi value at each point
  
  assert(xpt->cx1 < xpt->cx2);
  assert(xpt->cy1 < xpt->cy2);
  assert(xpt->cx1 >= 0);
  assert(xpt->cy1 >= 0);
  assert(xpt->cx2 < equilib->nw);
  assert(xpt->cy2 < equilib->nh);

  x = xpt->cx1;
  y = xpt->cy1;
  d = 1;
  n = 0;

  while (n < 6)
  {
    lvl = equilib->psi[x][y];
    
    ox = x;
    oy = y;

    if (y == xpt->cy1)
      x == xpt->cx2 ? y++ : x++;
    else if (x == xpt->cx2)
      y == xpt->cy2 ? x-- : y++;
    else if (y == xpt->cy2)
      x == xpt->cx1 ? y-- : x--;
    else if (x == xpt->cx1)
      y == xpt->cy1 ? x++ : y--;
    else
      assert(0);

    lvl1 = equilib->psi[x][y];

    if (d == 1 && lvl1 < lvl)
    {
      p[n].x = ox;
      p[n].y = oy;
      p[n].t = d;
      p[n].lvl = lvl;
      n++;
      d = -1;
    }
    else if (d == -1 && lvl1 > lvl)
    {
      p[n].x = ox;
      p[n].y = oy;
      p[n].t = d;
      p[n].lvl = lvl;
      n++;
      d = 1;
    }

    if (x == xpt->cx1 && y == xpt->cy1 && n == 0)
      break;
  }
// no extrema, not the xpoint rectangular.

  if (!n)
  {
    return -1; 
  }

  assert(n == 6);

// Check for exactly 4 extrema */
  for (i = 2; i < n - 1; i++)
  {
    if (p[i].x == p[1].x && p[i].y == p[1].y)
    {
      return -1; /* Too few */
    }
  }

  if (p[n - 1].x != p[1].x || p[n - 1].y != p[1].y)
  {
    return -1; /* Too many */
  }

  /* Make sure a minimum is 1st */

  if (p[1].t != -1)
    for (i = 1; i < 4; i++)
    {
      p[0] = p[i];
      p[i] = p[i + 1];
      p[i + 1] = p[0];
    }

  /* Shift values to the left */
  for (i = 0; i < 4; i++)
    p[i] = p[i + 1];

  n = 4;

  /* Check for each minimum < each maximum */
  if (p[0].lvl >= p[1].lvl || p[0].lvl >= p[3].lvl)
    return -1;
  if (p[2].lvl >= p[1].lvl || p[2].lvl >= p[3].lvl)
    return -1;

  /* Check for an 'X'-intersection */
  for (i = 0; i < 4; i++)
  {
    if (check_xpt_levels(equilib, xpt->cx1, xpt->cy1, xpt->cx2, xpt->cy2, p[i].x, p[i].y, p[i].t))
      return -1;
  }

  xpt->lvlMin = max(p[0].lvl, p[2].lvl);
  // printf("p[0].lvl: %.8lf\n", p[0].lvl);
  // printf("p[2].lvl: %.8lf\n", p[2].lvl);
  // printf("xpt->lvlMin: %.8lf\n", xpt->lvlMin);
  // printf("\n");

  xpt->lvlMax = min(p[1].lvl, p[3].lvl);
  // printf("p[1].lvl: %.8lf\n", p[1].lvl);
  // printf("p[3].lvl: %.8lf\n", p[3].lvl);
  // printf("xpt->lvlMmax: %.8lf\n", xpt->lvlMax);
  // printf("\n");

  for (i = 0; i < 4; i++)
    xpt->minMax[i] = p[i];

  return 0;
}

static int check_xpt_est_range(const Equilibrium *equilib, int xpoint_number, double **est_xpoint_pos, const XPointInfo xpt)
{
  if (!equilib || !est_xpoint_pos || !xpt) 
  {
    printf("Invalid input!\n");
    return -1;
  }
  double r,z;
  r = (equilib->r[xpt->cx1] + equilib->r[xpt->cx2])/2.0;
  z = (equilib->z[xpt->cy1] + equilib->z[xpt->cy2])/2.0;

  for (int i=0; i<xpoint_number; i++)
  {
    //printf("debug: r-est_xpoint_pos[i][0] %.12f\n", r-est_xpoint_pos[i][0]);
    //printf("debug: z-est_xpoint_pos[i][1] %.12f\n", z-est_xpoint_pos[i][1]);
    // printf("debug: r %.12f\n", r);
    // printf("debug: z %.12f\n", z);
    if (fabs(r-est_xpoint_pos[i][0])<DELTA && fabs(z-est_xpoint_pos[i][1])<DELTA)
    {
        return 0;
    }
  }
  return -1;
}

// refer DivGeo source code [static int CheckXPointLevels(Equil eq,int cx1,int cy1,int cx2,int cy2,int x0,int y0,int bMinMax)]
static int check_xpt_levels(const Equilibrium *equilib,int cx1,int cy1,int cx2,int cy2,int x0,int y0,int bMinMax)
{
  int x, y;
  int nx = 0, ny = 0;
  int d, nd;

  if ((x0 == cx1 || x0 == cx2) && (y0 == cy1 || y0 == cy2)) {
/*    AddSource(a,eq->x[x0],eq->y[y0]); */
    return 0;
  }
  x = x0;
  y = y0;
  
  if (x==cx1) {d=1;x++;} 
  else if (x==cx2) {d=3;x--;} 
  else if (y==cy1) {d=2;y++;} 
  else if (y==cy2) {d=0;y--;} 
  else assert(0);

  if((equilib->psi[x0][y0]-equilib->psi[x][y]) * bMinMax > 0 )
    return 0;

  do {
    for (nd=d-1;nd<=d+2;nd++) {
      nx=x;
      ny=y;
      switch(nd & 3) {
        case 0:ny=y-1;break;
        case 1:nx=x+1;break;
        case 2:ny=y+1;break;
        case 3:nx=x-1;break;
        default:assert(0);
      }

      if ((equilib->psi[x0][y0]-equilib->psi[x][y]) * bMinMax <= 0) 
      {
        d=nd;
        break;
      }
    }

    x=nx;
    y=ny;
/*printf("%d %d    ",x,y);*/
  }while (nx!=cx1 && nx!=cx2 && ny!=cy1 && ny!=cy2);
/*  if (nx==x0 && ny==y0) puts("Yau!"),AddSource(a,eq->x[x0],eq->y[y0]);*/
/*  else puts("AAA"); */
  return nx!=x0 || ny!=y0;
}


//Calculate the coordinates and psi values of the X-points 
static void calculate_xpt_level(Equilibrium *equilib, int xpoint_number, double **est_xpoint_pos,
                                interpl_1D_f interpl_1D_f, interpl_2D_f interpl_2D_f, _XPointInfo *xpt_array)
{
//currently, assume in the range [cx1:cx2,cy1:cy2] the psi are calculated by bilinear interpolation.

// for(int i=0; i<xpoint_number; i++)
// {
//   int cx1 = xpt_array[i].cx1;
//   int cy1 = xpt_array[i].cy1;
//   int cx2 = xpt_array[i].cx2;
//   int cy2 = xpt_array[i].cy2;
//   double a=(equilib->psi[cx2][cy1] - equilib->psi[cx1][cy1])/(equilib->r[cx2]-equilib->r[cx1]);
//   double b=(equilib->psi[cx1][cy2] - equilib->psi[cx1][cy1])/(equilib->z[cy2]-equilib->z[cy1]);
//   double c = (equilib->psi[cx2][cy2] + equilib->psi[cx1][cy1] - equilib->psi[cx2][cy1] - equilib->psi[cx1][cy2])
//             /((equilib->r[cx2]-equilib->r[cx1])*(equilib->z[cy2]-equilib->z[cy1]));

//   xpt_array[i].centerX = (-b)/c + equilib->r[cx1];
//   xpt_array[i].centerY = (-a)/c + equilib->z[cy1];
//   xpt_array[i].level = equilib->psi[cx1][cy1] - (a*b)/c;
//   printf("R: %.12f, Z: %.12f, psi: %.12f\n", xpt_array[i].centerX, xpt_array[i].centerY, xpt_array[i].level);
//   printf("crx1: %i, crx2: %i, cry1: %i, cry2: %i\n", cx1, cx2, cy1, cy2);
//   printf("x1: %.12f y1: %.12f\n", equilib->r[cx1], equilib->z[cy1]);
//   printf("lvlmin: %.12f, lvlmax: %.12f\n", xpt_array[i].lvlMin, xpt_array[i].lvlMax);
// }

/******************************************************************************
* Following is another way based on secant method. There are problems in that 
* method. In the future, it will be used.
*******************************************************************************/
  //temperoray magneticfield used for calculate X-point where the magnetic field is zero.

  MagFieldTorSys magfield;
  init_mag_field_torsys(&magfield);
  char* method = "central_4th";
  calc_mag_field_torsys(equilib, &magfield, method);
  //write_mag_field_torsys(&test_magfield);

  // //temperoray magneticfield used for calculate X-point where the magnetic field is zero.
  // double ***grad_psi = allocate_3d_array(equilib->nw, equilib->nh, 2);
  // central_2nd_2d_diff(equilib->nw, equilib->r, equilib->nh, equilib->z, equilib->psi, grad_psi);


  for(int i=0; i<xpoint_number; i++)
  {
    int cx1 = xpt_array[i].cx1;
    int cy1 = xpt_array[i].cy1;
    int cx2 = xpt_array[i].cx2;
    int cy2 = xpt_array[i].cy2;
    
    double x_start = equilib->r[cx1];
    double y_start = equilib->z[cy1];
    double x = x_start;
    double y = y_start;
    double dBrdx = magfield.dBrzdx[cx1][cy1][0];
    double dBrdy = magfield.dBrzdy[cx1][cy1][0];
    double dBzdx, dBzdy;
    double L2norm;
    double Br,Bz;
    interpl_2D_f(x,y, equilib->nw, equilib->r, equilib->nh, equilib->z,
                    magfield.Brz, &Br, &Bz, NULL, NULL, NULL);
    
    L2norm = sqrt(pow(Br,2.0) + pow(Bz,2.0));
    printf("x: %.12f y: %.12f L2norm: %.12f\n", x, y, L2norm);
    printf("Br: %.12f Bz: %.12f\n", Br, Bz);
    if(L2norm < TOLERANCE)
    {
      xpt_array[i].centerX = x_start;
      xpt_array[i].centerY = y_start;
      continue;
    }

    //printf("x: %.12f y: %.12f\n", x, y);
    for (int iter=0; iter<MAXITER; iter++)
    {
      interpl_2D_f(x,y, equilib->nw, equilib->r, equilib->nh, equilib->z,
                  magfield.dBrzdx, &dBrdx, &dBzdx, NULL, NULL, NULL);

      interpl_2D_f(x,y, equilib->nw, equilib->r, equilib->nh, equilib->z,
                  magfield.dBrzdy, &dBrdy, &dBzdy, NULL, NULL, NULL);



      //we always start from the left corner. The Xpoint is in the center.
      //In order to make sure the search in the correct direction,
      // we use x+LEARNINGRATE*(fabs(dBrdx)+fabs(dBzdx))/2.0 instead of x-LEARNINGRATE*dBrdx
      // also y=y+LEARNINGRATE*(fabs(dBrdx)+fabs(dBzdx))/2.0 instead of y-LEARNINGRATE*dBrdy
      // At the meanwhile, the the effec of dBzdx and dBzdy is considered.


      x=x+LEARNINGRATE*(fabs(dBrdx)+fabs(dBzdx))/2.0;
      y=y+LEARNINGRATE*(fabs(dBrdy)+fabs(dBzdy))/2.0;


      interpl_2D_f(x,y, equilib->nw, equilib->r, equilib->nh, equilib->z,
                  magfield.Brz, &Br, &Bz, NULL, NULL, NULL);
    
      L2norm = sqrt(pow(Br,2.0) + pow(Bz,2.0));
      //printf("iter: %d\n", iter);
      //printf("x: %.12f y: %.12f L2norm: %.12f\n", x, y, L2norm);
      //printf("Br: %.12f Bz: %.12f\n", Br, Bz);

      //printf("dfdx: %.12f dfdy: %.12f\n", dBrdx, dBrdy);

      if(L2norm < TOLERANCE)
      {
        xpt_array[i].centerX = x;
        xpt_array[i].centerY = y;
        printf("iter: %d\n", iter);
        printf("x: %.12f y: %.12f L2norm: %.12f\n", x, y, L2norm);
        printf("Br: %.12f Bz: %.12f\n", Br, Bz);
        break;
      }
      if(iter == MAXITER-1)
      {
        printf("Don't find the X-point center. Please check the function: calculate_xpt_level.\n");
        exit(EXIT_FAILURE);
      }
    }

    //double check the xpoint
    if (xpt_array[i].centerX >= equilib->r[cx1] &&
        xpt_array[i].centerX <= equilib->r[cx2] &&
        xpt_array[i].centerY >= equilib->z[cy1] &&
        xpt_array[i].centerY <= equilib->z[cy2])
    {
      printf("The %d X-point is in the range\n", i);
      printf("The X-Point: %.12f %.12f\n", xpt_array[i].centerX, xpt_array[i].centerY);
    }
    else
    {
      printf("The %d X-point is in out of the range\n", i);
      printf("Please check the function: calculate_xpt_level\n");
      exit(1);
    }
  }

  for(int i=0; i<xpoint_number; i++)
  {
    interpl_1D_f(xpt_array[i].centerX, xpt_array[i].centerY, 
                equilib->nw,equilib->r, equilib->nh, equilib->z, equilib->psi,
                &(xpt_array[i].level), NULL, NULL, NULL);
    printf("The psi of the %d X-Point: %.12f \n", i, xpt_array[i].level);
  }
  free_mag_field_torsys(&magfield);
}

void test_find_xpoint()
{
  printf("*******************************************\n");
  printf("Begint to test find_xpoint\n"); 

  InputPara w3_input;
  init_inputpara(&w3_input);
  print_inputpara(&w3_input);
  Equilibrium dtt_example;
  
  init_equilibrium(&dtt_example);

  read_equilib_geqdsk(&dtt_example,w3_input.equilibrium_file);
  print_equilibrium(&dtt_example);
  
  int xpt_n = 2;
  double **est_xpt = allocate_2d_array(xpt_n,2);
  est_xpt[0][0] = 1.85;
  est_xpt[0][1] = -1.16;

  est_xpt[1][0] = 1.58;
  est_xpt[1][1] = 1.61;

  interpl_1D_f interpl_1D_f = bilenar_1d;
  interpl_2D_f interpl_2D_f = bilenar_2d;

  _XPointInfo xpt_array[2];

  find_xpoint(&dtt_example, xpt_n, est_xpt, interpl_1D_f, interpl_2D_f, xpt_array);

  free_2d_array(est_xpt);
  free_equilibrium(&dtt_example);  

}