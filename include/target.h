#ifndef TARGET_H
#define TARGET_H

#include "datastructure.h"
#include "divgeo.h"

typedef struct 
{
  char name[32];
  int n; //number of points
  DLListNode* head;
} TargetCurve;

TargetCurve* create_target_curve();
//n indicate which target curve in trg
TargetCurve* create_target_curve_from_dgtrg(DivGeoTrg* trg, int n);

int add_point_target_curve(TargetCurve* curve, double r, double z);

int change_name_target_curve(TargetCurve* curve, const char* name);

void free_target_curve(TargetCurve* curve);
//create a new target curve which has the opposite sequence.

//pass the adress of a TargetCurve*, the old curve is free and a new 
//curve is created and bound to the old pointer
TargetCurve* reverse_target_curve(TargetCurve** curve);

void printf_target_curve(TargetCurve* curve);





#endif