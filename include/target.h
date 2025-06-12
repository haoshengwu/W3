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

int add_point_target_curve(TargetDLListCurve* tgt_cur, double r, double z);

int change_name_target_curve(TargetDLListCurve* tgt_cur, const char* name);

void free_target_curve(TargetDLListCurve* tgt_cur);

//split the DDList in a target curve and the new DDList will be bound to new_tgt_cur
void split_intersections_target_curve(TargetDLListCurve* tgt_cur, 
                                      double r, double z, 
                                      TargetDLListCurve* new_tgt_cur);


//update the number of target curve
void update_number_target_curve(TargetDLListCurve* tgt_cur);

//reverse the DLList in target_cureve
void reverse_DDList_in_target_curve(TargetDLListCurve* tgt_cur);

void printf_target_curve(TargetDLListCurve* tgt_cur);


//the sep line which intersect with target curve is 1.
//gradpsi line is between sepline[i] and sepline[i+1]
//
//                 \    :     /
//                  \   :    /
//            ...... X-point......gradpsi line
//                  /   :   \
//                 /    :    \ sep line

void sort_sep_gradpsiline_by_targetcurve(TargetDLListCurve* tgt_cur, SeparatrixStr* sep, GradPsiLineStr* gradpsilines);





#endif