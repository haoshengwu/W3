#ifndef XPOINT_H
#define XPOINT_H
#include "equilibrium.h"
#include "mathbase.h"



// Some structs are from DG, the objects with common out are not used.

//Store the info of the extremums (maxmium and minium) for a Xpoint.
//x and y are the point number, t indicates is maximum or minium.
//lvl is the psi value;
struct XPointExtremum {
  int x,y,t;
  double lvl;
};

typedef struct {
  int type,locks;
  //App app;
  int cx1,cy1,cx2,cy2;
  double lvlMin,lvlMax,level,centerX,centerY;
  struct XPointExtremum minMax[4];
  //Group segs;
  //Group gradients[4];
  //int id;
} _XPointInfo;

typedef _XPointInfo* XPointInfo;


void find_xpoint(Equilibrium *equ, int xpoint_number, double **est_xpoint_pos,
                 interpl_2D_1f interpl_2D_1f, interpl_2D_2f interpl_2D_2f, _XPointInfo *xpoint);


void test_find_xpoint();
#endif