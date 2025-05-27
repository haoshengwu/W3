#include "magneticsurface.h"



#define CSF_YP    0x1
#define CSF_YM    0x2
#define CSF_XM    0x4
#define CSF_XP    0x8

#define CS_YP    0
#define CS_XP    1
#define CS_YM    2
#define CS_XM    3


#define EqXEx(eq,cx,_nw)\
  ((eq)->r[0]+((eq)->r[(eq)->nw-1]-(eq)->r[0])*(cx)/((_nw)-1))

#define EqYEx(eq,cy,_nh)\
  ((eq)->z[0]+((eq)->z[(eq)->nh-1]-(eq)->z[0])*(cy)/((_nh)-1))

//This function calcuate the line for one cell, accroding DG [CalcSurfData in equil.c]
//sx, and sy are the number of points, not the total cell number.
void calc_surf_data(const Equilibrium *equlib,int cx,int cy,double level,struct _SurfCell* sc,int sx,int sy)
{
  double a1, a2;
  int d;
  sc->n = sc->f = 0;
  if (cx<0 || cy<0 || cx>=sx-1 || cy>=sy-1) 
  {
    printf("arriave at the boundary cx or cy is out of range\n");
    return;
  }
  // currently use EqCorrCell, should be consist with DG use EqCorrCellEx
  a1 = EqCorrCell(equlib, cx,cy,level);
  a2 = EqCorrCell(equlib, cx+1,cy,level);
  if (inrange_s(level, a1,a2))
  {
    d = CS_YM;
    sc->x[d] = EqXEx(equlib, cx, sx) + (EqXEx(equlib,cx+1,sx)-EqXEx(equlib,cx,sx))*(level-a1)/(a2-a1);
    sc->y[d] = EqYEx(equlib, cy, sy);
    sc->f |= CSF_YM;
    sc->d[sc->n] = d;
    sc->n++;
  }

  a2=EqCorrCell(equlib,cx,cy+1,level);
  if (inrange_s(level,a1,a2))
  {
    d=CS_XM;
    sc->x[d]=EqXEx(equlib,cx,sx);
    sc->y[d]=EqYEx(equlib,cy,sy)+
        (EqYEx(equlib,cy+1,sy)-EqYEx(equlib,cy,sy))*(level-a1)/(a2-a1);
    sc->f |= CSF_XM;
    sc->d[sc->n++]=d;
  }
  a1=EqCorrCell(equlib,cx,cy+1,level);
  a2=EqCorrCell(equlib,cx+1,cy+1,level);
  if (inrange_s(level,a1,a2)) 
  {
    d=CS_YP;
    sc->x[d]=EqXEx(equlib,cx,sx)+
        (EqXEx(equlib,cx+1,sx)-EqXEx(equlib,cx,sx))*(level-a1)/(a2-a1);
    sc->y[d]=EqYEx(equlib,cy+1,sy);
    sc->f |= CSF_YP;
    sc->d[sc->n++]=d;
  }
  a1=EqCorrCell(equlib,cx+1,cy,level);
  if (inrange_s(level,a1,a2)) 
  {
    d=CS_XP;
    sc->x[d]=EqXEx(equlib,cx+1,sx);
    sc->y[d]=EqYEx(equlib,cy,sy)+
        (EqYEx(equlib,cy+1,sy)-EqYEx(equlib,cy,sy))*(level-a1)/(a2-a1);
    sc->f |= CSF_XP;
    sc->d[sc->n++]=d;
  }
  return;
};


// This function is according to [CalcSurfaceLine] in DG. The input are: the starting cell number [cx,cy]
// The level is the prescribe psi value. This function will trace the line which have the same psi [level]
// In the future, it will be replaced by a file line tracer like FLARE and INGREID.
/* Calculate group of points representing a surface
   Return: -1 = invalid origin
            0 = ok, surface not closed
            1 = ok, surface closed
*/
// All the coordinate r,z of the points are stored in a double-lined list 'ptr_line_list'
// nw and nh are total number of points, not the cell number.
int calc_surface_line(const Equilibrium *equlib,int cx,int cy,double level,int nw,int nh, DLListNode** ptr_line_list)
{
  struct _SurfCell sc;
  int oCx,oCy,d,closed=0;
  oCx = cx;
  oCy = cy;
  DLListNode *endnode = NULL;


  if ((*ptr_line_list) != NULL)
  {
    printf("!!!the DDListNod in not empty, are you sure????\n");
  }

  calc_surf_data(equlib, cx, cy, level, &sc, nw, nh);
  if (sc.n!=2 && sc.n!=4) 
  {
    printf("in this cell cx:%d, cy:%d, there is no prescribed value %lf", cx,cy,level);
    return -1;
  }

  d = sc.d[0];

  do {
    // out put the intersection points, later will be change
    insert_DLList_at_head(ptr_line_list, sc.x[d], sc.y[d]);

    // printf("%lf %lf\n", sc.x[d], sc.y[d]);
    if (d==CS_YM) cy--;
    if (d==CS_YP) cy++;
    if (d==CS_XM) cx--;
    if (d==CS_XP) cx++;
    calc_surf_data(equlib, cx, cy, level, &sc, nw, nh);
    
    if (sc.n==2) 
    {
      // the boundary in previous cell is opposite of the current cell
      if ( (d^2) == sc.d[0]) 
        d = sc.d[1];
      else 
        d = sc.d[0];
    }
    // back to the origin point, so it is closed
    if (cx==oCx && cy==oCy) 
    {
      printf("back to the origin, cx = %d, cy = %d", cx, cy);
      closed = 1;
      return closed;
    }
  } while (sc.n == 2 || sc.n == 4);

  // go back to the beginning and trace in another direction
  cx = oCx;
  cy = oCy;
  calc_surf_data(equlib, cx, cy, level, &sc, nw, nh);
  
  endnode = get_DLList_endnode(*ptr_line_list);
  //printf("debug: first endnode is %p, r: %lf, z: %lf\n", (void*)endnode, endnode->r,endnode->z);

  if (sc.n == 2)
    d = sc.d[1];
  else
    d = sc.d[0]^2;
  do
  {
    insert_DLList_at_end(&endnode, sc.x[d], sc.y[d]);
    //printf("debug: the endnode is %p, r: %lf, z: %lf\n", (void*)endnode, endnode->r,endnode->z);

    // out put the intersection points, later will be change
    // printf("%lf %lf\n", sc.x[d], sc.y[d]);
    if (d==CS_YM) cy--;
    if (d==CS_YP) cy++;
    if (d==CS_XM) cx--;
    if (d==CS_XP) cx++;
    calc_surf_data(equlib, cx, cy, level, &sc, nw, nh);

    if (sc.n==2) 
    {
      // the boundary in previous cell is opposite of the current cell
      if ( (d^2) == sc.d[0]) 
        d = sc.d[1];
      else 
        d = sc.d[0];
    }
  } while (sc.n == 2 || sc.n == 4);

  return closed;

};

//void cal_magnetic_surface(const Equilibrium *equlib, const double psi_level, MagneticSurfaceSegment segent)
//{
//  return;
//};


// Refer to DivGeo [CalcSeparatrixLine] in [xpoint.c]
// idx indicates which point is the start.
// A X-point have four segments, the idx indicate which is the one used. 
void cal_separatrix_line(const Equilibrium *equlib, const XPointTest xpt, int idx, DLListNode** ptr_line_list)
{
  
  int n,x,y,ox,oy;
  // struct _SurfCell sc;
  // XY xy,xy1,xy0;

// check the xpt positions
  assert(xpt->cx1 > 0);
  assert(xpt->cy1 > 0);
  assert(xpt->cx2 < equlib->nw-1);
  assert(xpt->cy2 < equlib->nh-1);
  assert(xpt->cx1 < xpt->cx2);
  assert(xpt->cy1 < xpt->cy2);

  n=0;

  x=xpt->cx1;
  y=xpt->cy1;

  // 
  while(1) {
    ox=x;
    oy=y;
    if (y==xpt->cy1) x==xpt->cx2 ? y++ : x++; 
    else if (x==xpt->cx2) y==xpt->cy2 ? x-- : y++; 
    else if (y==xpt->cy2) x==xpt->cx1 ? y-- : x--; 
    else if (x==xpt->cx1) y==xpt->cy1 ? x++ : y--; 
    else assert(0);

    if (inrange_s(xpt->level,EqCorrCell(equlib,ox,oy,xpt->level),EqCorrCell(equlib,x,y,xpt->level)))
      if (n++==idx) break;

    if (x==xpt->cx1 && y==xpt->cy1) 
    {/* puts("1"); */ printf("please check x-point");}
  }

  // here x and y is the cell number, not the point number.
  //
  //      |   cell  | 
  //---------------------
  //      |         |
  // cell | X-point |  cell
  //      |         |
  //--------------------
  //      |   cell  |
  //
  if (x>ox) swap(x,ox);
  if (y>oy) swap(y,oy);
  
  // Again, adjust the cell number in different situation
  // tracing not from the X-point rectangular, but from the neighbour cells.

  if (x==ox) 
  {
      if (x==xpt->cx1) 
          x--;
  } else if (y==oy) 
  {
    if (y==xpt->cy1) 
        y--;
  } 
  else {
    assert(0);
  }
  
    int i;
    printf("tracing the segment %d of separatrix:\n", idx); 
    i = calc_surface_line(equlib, x, y, xpt->level, equlib->nw, equlib->nh, ptr_line_list);
    printf("the situation of line is: %d\n",i);

}
