#ifndef MAGNETICSURFACE_H
#define MAGNETICSURFACE_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "equilibrium.h"
#include "calc.h"
#include "datastructure.h"

#define CSF_YP    0x1
#define CSF_YM    0x2
#define CSF_XM    0x4
#define CSF_XP    0x8

#define CS_YP    0
#define CS_XP    1
#define CS_YM    2
#define CS_XM    3



#define POINTS_LIMIT 2000 // the maxium number of points in a magnetic surface to store the positions.

struct _MagneticSurfaceSegment
{
    double psi_level;
    double r[POINTS_LIMIT];
    double z[POINTS_LIMIT];
    int close;
    int intersetion_num;
};

struct _SurfCell {
  // r and z coordinates of the intersection points for the four boundaries of a cell.
  double x[4],y[4]; 
  // n is the total intersection number, 
  // f indicates the situations, there is total 16 situations. [need double check]
  int n,f;
  // indicates the boundary numberi. e.g. d[0] = 2 meanse the first intersection point is at the 'CS_YM = 2' boundary. 
  int d[4];
};

struct _XY {
  double x,y;
};

typedef struct _XY* XY;

typedef struct _MagneticSurfaceSegment* MagneticSurfaceSegment;

void calc_surf_data(const Equilibrium *equlib,int cx,int cy,double level,struct _SurfCell* sc,int sx,int sy);

int calc_surface_line(const Equilibrium *equlib,int cx,int cy,double level,int nw,int nh, DLListNode** ptr_line_list);

void cal_magnetic_surface(const Equilibrium *equlib, const double psi_level, MagneticSurfaceSegment segment);

void cal_separatrix_line(const Equilibrium *equlib, const XPointTest xpt, int index, DLListNode** ptr_line_list );

void write_magnetic_surface(const MagneticSurfaceSegment segment, char *filename);



#endif
