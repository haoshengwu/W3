#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "equilibrium.h"
#include "calc.h"

// private functions for find Xpoint
static int check_xpt_rectangular(const Equilibrium *equilib, const XPointTest xpt);
static int check_point_in_rectangular(const double pos[], const Equilibrium *equilib, int cx1, int cx2, int cy1, int cy2);
static int check_xpt_levels(const Equilibrium *equilib,int cx1,int cy1,int cx2,int cy2,int x0,int y0,int bMinMax);
static void calculate_xpt_level(XPointTest xpt);
static int calculate_xpt_center(Equilibrium *equilib, XPointTest xpt);


//refer to DivGep [void CalcEquilValues(Equil eq)]
void calculate_equi_values(Equilibrium *equilib)
{
  int i,j;
  equilib->minVal = equilib->psi[0][0];
  equilib->maxVal = equilib->psi[0][0];
  for(i=0; i < equilib->nw; i++)
    for (j= 0; i < equilib->nh; j++)
    {
      equilib->minVal = min(equilib->minVal, equilib->psi[i][j]);
      equilib->minVal = max(equilib->maxVal,equilib->psi[i][j]);
    }
}

double EqCorrCell(const Equilibrium *equilib,int cx,int cy,double level) 
{
  double a;
  a = equilib->psi[cx][cy];
  return a==level ? a + (equilib->maxVal - equilib->minVal)*1e-7 : a;
}

void init_equilibrium(Equilibrium *equilib)
{
  equilib->nw = 0;
  equilib->nh = 0;
  equilib->simag = 0;
  equilib->sibry = 0;
  equilib->r = NULL;
  equilib->z = NULL;
  equilib->psi = NULL;
  
  equilib->Xpoint_num = 0;
  equilib->Xpoint_pos[0] = 0;
  equilib->Xpoint_pos[1] = 0;
  equilib->Opoint_num = 0;
  equilib->Opoint_pos[0] = 0;
  equilib->Opoint_pos[1] = 0;

  equilib->minVal = 0;
  equilib->maxVal = 0;
}

void print_equilibrium(const Equilibrium *equilib)
{
  printf("Equilibrium size:\n");
  printf("nw: %i, nh: %i\n", equilib->nw, equilib->nh);
  printf("poloidal flux at magnetic axis is %.8lf Weber/rad\n", equilib->simag);
  printf("poloidal flux at the plasma boundary is %.8lf Weber/rad\n", equilib->sibry);

  // print part of equilibrirum value to check
  // for (int i = 0; i < 5; i++)
  // {
  //   printf("psi i = %i, j = 1: %.8lf\n", i + 1, equilib->psi[i][0]);
  // }
  // for (int i = 0; i < 5; i++)
  // {
  //   printf("psi i = %i, j = %i: %.8lf\n", i + 1, equilib->nh, equilib->psi[i][equilib->nh - 1]);
  // }
}

//  Currently only store necessary paramters/
void read_equilib_geqdsk(Equilibrium *equilib, const char *geqdsk_file)
{

  FILE *file = fopen(geqdsk_file, "r");

  if (file == NULL)
  {
    fprintf(stderr, "Error opening file %s\n", geqdsk_file);
    exit(1); // Failed to read geqdsk equilibrium;
  }

  //  the length of title is 48 + 1 for /0;
  //
  char title[49];
  int idum, nw, nh;
  double rdim, zdim, rcentr, rleft, zmid;
  double rmaxis, zmaxis, simag, sibry, bcentr;
  double current, xdum;
  double tmp;

  fscanf(file, "%48s %i %i %i", title, &idum, &nw, &nh);
  fscanf(file, "%lf %lf %lf %lf %lf", &rdim, &zdim, &rcentr, &rleft, &zmid);
  fscanf(file, "%lf %lf %lf %lf %lf", &rmaxis, &zmaxis, &simag, &sibry, &bcentr);
  fscanf(file, "%lf %lf %lf %lf %lf", &current, &simag, &xdum, &rmaxis, &xdum);
  fscanf(file, "%lf %lf %lf %lf %lf", &zmaxis, &xdum, &sibry, &xdum, &xdum);

  equilib->nw = nw;
  equilib->nh = nh;
  equilib->simag = simag;
  equilib->sibry = sibry;

  // skip fpol, pres, ffprim, pprime, total number is 4 * nw, will be updated.
  for (int i = 0; i < 4 * equilib->nw; i++)
  {
    fscanf(file, "%lf", &tmp);
  }

  //  Allocate dynamic memmory and passing value for equilib -> psi
  equilib->psi = (double **)malloc(equilib->nw * sizeof(double *));

  for (int i = 0; i < equilib->nw; i++)
  {
    equilib->psi[i] = (double *)malloc(equilib->nh * sizeof(double));
  }
  
 
  for (int j = 0; j < equilib->nh; j++)
  {
    for (int i = 0; i < equilib->nw; i++)
    {
      fscanf(file, "%lf", &equilib->psi[i][j]);
    }
  }

  printf("finish reading geqdsk\n");
  // Allocate memmory and calculate r and z values.

  equilib->r = (double *)malloc(equilib->nw * sizeof(double));
  equilib->z = (double *)malloc(equilib->nh * sizeof(double));

  double delta_r = rdim / (equilib->nw - 1);
  for (int i = 0; i < equilib->nw; i++)
  {
    equilib->r[i] = i * delta_r + rleft;
  }

  double delta_z = zdim / (equilib->nh - 1);
  double zleft = (2 * zmid - zdim) /2 ;
  for (int j = 0; j < equilib->nh; j++)
  {
    equilib->z[j] = zleft + j * delta_z;
  }
  printf("r, z and psi is ok\n");

}

XPointTest find_Xpoint(Equilibrium *equilib, const double *est_pos)
{
  // est_pos is the estimated postion est_pos[0] is R, est_pos[1] is Z. n is the size of est_pos
  // this algorithm is from DivGeo FindXPointRects and FindXPointCenter xpoint.h and xpoint.c files.
  int i, j, s;

  XPointTest xp, xpC; // xpC is for every potentional xpoint. xp is the one store the X-point.

  xpC = malloc(sizeof(*xpC));
  xp = malloc(sizeof(*xp));
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

        // printf("Find a potential equ rectangular contain a Xpoint:\n");
        // printf("crx1: %i, crx2: %i, cry1: %i, cry2: %i\n", xpC->cx1, xpC->cx2, xpC->cy1, xpC->cy2);
        // printf("R: %.8lf Z: %.8lf\n", xpC->centerX, xp->centerY);

        // store the value in xp
        if (!check_point_in_rectangular(est_pos, equilib, xpC->cx1,xpC->cx2,xpC->cy1,xpC->cy2))
        {
          xp->cx1 = xpC->cx1;
          xp->cx2 = xpC->cx2;
          xp->cy1 = xpC->cy1;
          xp->cy2 = xpC->cy2;
          xp->lvlMin = xpC->lvlMin;
          xp->lvlMax = xpC->lvlMax;
        }
      }
    }
  }

  calculate_xpt_level(xp);
  // printf("calculate_xpt_level: %.8lf\n",xp->level);
  
 if(calculate_xpt_center(equilib,xp))
  {
    printf("NOT find the X point, PLEASE check");
  };

  printf("crx1: %i, crx2: %i, cry1: %i, cry2: %i\n", xp->cx1, xp->cx2, xp->cy1, xp->cy2);
  printf("The X-point is R: %.8lf, Z: %.8lf\n", xp->centerX, xp->centerY);
  printf("r range: %.8lf, %.8lf\n", equilib->r[xp->cx1], equilib->r[xp->cx2]);
  printf("z range: %.8lf, %.8lf\n", equilib->z[xp->cy1], equilib->z[xp->cy2]);
  printf("value: %.8lf\n", xp->level);

  // store the value in equlib structure
  equilib->Xpoint_num = 1;
  equilib->Xpoint_pos[0] = xp->centerX;
  equilib->Xpoint_pos[1] = xp->centerY;
  
  return xp; // !!!DO NOT free xp and xpC
  free(xpC);
  // free(xp);
}

// x is the r coordinate, y is the z coordinate, we assume x,y is not any apoint of equilibrium
double get_psi_from_rz(const Equilibrium *equilib, double x, double y)
{
  // value is the psi at the point (x,y)
  double value;
  if( x > equilib->r[equilib->nw - 1] || x < equilib->r[0])
  {
    printf("the position of x: %lf is out of r range\n", x);
    return NAN;
  }

  if( y > equilib->z[equilib->nh - 1] || y < equilib->z[0])
  {
    printf("the position of y: %lf is out of z range\n", y);
    return NAN;
  }
  
  // this is the cell number, not the point. the coresponding four corners are: (cx,cy), (cx+1,cy), (cx, cy+1), (cx+1, cy+1) 
  int cx, cy;
  
  for (int j = 0; j < equilib->nh; j++)
  {
    if ( y > equilib->z[j] && y < equilib->z[j+1])
    {
      for (int i = 0; i < equilib->nw; i++)
        {
          if ( x > equilib->r[i] &&  x < equilib->r[i+1])
          {
            cx = i;
            cy = j;
            printf("the x: %lf, y: %lf in the cell cx: %d cy: %d\n",x,y,cx,cy);
            break;
          }
        }
    }
  }
  
  // Bi-linear interpolation for the psi valuex according Kotov's [Generation of orthogonal full-device grid]
  //psi(x,y) = psi[cx][cy] + a*(x-r[cx]) + b*(y-z[cy]) + c*(x-r[cx])*(y-z[cy])
  double a,b,c;
  double delta_r,delta_z;
  delta_r = x - equilib->r[cx];
  delta_z = y - equilib->z[cy];
  // printf("%lf %lf\n",equilib->psi[cx][cy], equilib->psi[cx+1][cy]);
  // printf("%lf %lf\n",delta_r, delta_z);

  a = ( equilib->psi[cx+1][cy] - equilib->psi[cx][cy] ) / (equilib->r[cx+1] - equilib->r[cx] );
  b = ( equilib->psi[cx][cy+1] - equilib->psi[cx][cy] ) / (equilib->z[cy+1] - equilib->z[cy] );
  c = ( equilib->psi[cx+1][cy+1] + equilib->psi[cx][cy] - equilib->psi[cx+1][cy] - equilib->psi[cx][cy+1])/ \
      (equilib->r[cx+1] - equilib->r[cx]) / (equilib->z[cy+1] - equilib->z[cy]);
  value = equilib->psi[cx][cy] + a * delta_r + b*delta_z + c * delta_r * delta_z;
  return value;
}

// refer to DivGeo source code [static int CheckXPointRect(Equil eq,XPointTest xpt)]
static int check_xpt_rectangular(const Equilibrium *equilib, const XPointTest xpt)
{
  struct _XPointMinMax p[6];
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

// check whether estimate point in the rectangluar
static int check_point_in_rectangular(const double pos[], const Equilibrium *equilib, int cx1, int cx2, int cy1, int cy2)
{
  if ((pos[0] > equilib->r[cx1] && pos[0] < equilib->r[cx2]) \
  && (pos[1] > equilib->z[cy1] && pos[1] < equilib->z[cy2]))
  {
    // printf("inside?");
    return 0; // 0 menas the point pos[] in the rectangular;
  }
  else
    // printf("%.8lf, %.8lf, %.8lf\n", pos[0], equilib->r[cx1], equilib->r[cx2]);
    // printf("%.8lf, %.8lf, %.8lf\n", pos[1], equilib->z[cy1], equilib->z[cy2]);
    // printf("\n");

    return 1; // 1 menas the point pos[] NOT in the rectangular;
    
}

//calculate psi value

static void calculate_xpt_level(XPointTest xpt)
{
  printf("calculate_xpt_level: xpt->lvlMin: %.8lf\n", xpt->lvlMin);
  printf("calculate_xpt_level: xpt->lvlMax: %.8lf\n", xpt->lvlMax);
  xpt->level = (xpt->lvlMin + xpt->lvlMax)/2;
}

// caculate the X point postion, refer to DivGep [FindXPointCenter(Equil eq,XPointTest xpt)]
static int calculate_xpt_center(Equilibrium *equilib, XPointTest xpt)
{
  double xs[4],ys[4],lvl,lvl1,r;
  int n=0,x,y,ox,oy;

  x=xpt->cx1;
  y=xpt->cy1;

    while(n<4) {
    lvl=EqCorrCell(equilib,x,y,xpt->level);

    ox=x;
    oy=y;

    if (y==xpt->cy1) x==xpt->cx2 ? y++ : x++;
    else if (x==xpt->cx2) y==xpt->cy2 ? x-- : y++; 
    else if (y==xpt->cy2) x==xpt->cx1 ? y-- : x--; 
    else if (x==xpt->cx1) y==xpt->cy1 ? x++ : y--; 
    else assert(0);

    lvl1=EqCorrCell(equilib,x,y,xpt->level);
    
    // printf("xpt->level: %.8lf\n", xpt->level);
    if ((lvl-xpt->level)*(lvl1-xpt->level)<0) {
      xs[n]=equilib->r[ox]+(equilib->r[x]-equilib->r[ox])*(xpt->level-lvl)/(lvl1-lvl);
      ys[n]=equilib->z[oy]+(equilib->z[y]-equilib->z[oy])*(xpt->level-lvl)/(lvl1-lvl);
/*      AddSource(a,xs[n],ys[n]); */
      n++;
    }

    if (x==xpt->cx1 && y==xpt->cy1 && n==0) break;
  }

  if (n<4) return -1;

  if (VIntersect(xs[0],ys[0],xs[2],ys[2],xs[1],ys[1],xs[3],ys[3],&r,NULL))
    return -1;

  xpt->centerX=xs[0]+(xs[2]-xs[0])*r;
  xpt->centerY=ys[0]+(ys[2]-ys[0])*r;

/*  AddSource(a,xpt->centerX,xpt->centerY); */
  return 0;
}

void free_equilibrium(Equilibrium *equilib)
{

  for (int i = 0; i < equilib->nw; i++)
  {
    free(equilib->psi[i]);
    equilib->psi[i] = NULL;
  }
  free(equilib->psi);
  equilib->psi = NULL;

  free(equilib->r);
  free(equilib->z);
  equilib->r = NULL;
  equilib->z = NULL;
  equilib->nw = 0;
  equilib->nh = 0;
}

// void find_Xpoint(const Equilibrium* equilib, const double estimate[2], double accurate[2]){
//
// };