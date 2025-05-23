#ifndef SPARATRIX_H
#define SPARATRIX_H

#include "datastructure.h"
#include "xpoint.h"
#include "equilibrium.h"
#include "ode.h"
#include "tracertest.h"


typedef struct Separatrix Separatrix;


typedef struct Separatrix
{
  double r,z; // X-point;
  double psi; //psi
  int index[4]; // determine the index of septrix
  int order; //1st, 2nd, 3rh, or fouth X-point
  DLListNode* line_list[4]={NULL}; //four lines
  SeparatrixOpt* opt;
}
Separatrix;

typedef struct SeparatrixOpt
{
  Separatrix* init();
  void (*free)(Separatrix* sep);
  void (*generate)(Separatrix* sep, _XPointInfo* xpoint, Equilibrium* equ, ode_function ode_func, ode_solver brk45_solver);
}
SeparatrixOpt;


#endif SPARATRIX_H