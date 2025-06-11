#ifndef TARGET_H
#define TARGET_H

#include "datastructure.h"
#include "divgeo.h"
#include "separatrix.h"
#include "gradpsi.h"

typedef struct 
{
  char name[32];
  int n; //number of points
  DLListNode* head;
} TargetDLListCurve;

TargetDLListCurve* create_target_curve();

//n indicate which target curve in trg
TargetDLListCurve* create_target_curve_from_dgtrg(DivGeoTrg* trg, int n);

int add_point_target_curve(TargetDLListCurve* curve, double r, double z);

int change_name_target_curve(TargetDLListCurve* curve, const char* name);

void free_target_curve(TargetDLListCurve* curve);
//create a new target curve which has the opposite sequence.

//pass the adress of a TargetDLListCurve*, the old curve is free and a new 
//curve is created and bound to the old pointer
TargetDLListCurve* reverse_target_curve(TargetDLListCurve** curve);

void printf_target_curve(TargetDLListCurve* curve);

void sort_sep_gradpsilin_by_targetcure(TargetDLListCurve** targetcurve, SeparatrixStr* sep, GradPsiLineStr* gradpsilines);





#endif