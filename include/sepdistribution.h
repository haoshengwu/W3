#ifndef SEPDISTRIBUTION_H
#define SEPDISTRIBUTION_H

#include "datastructure.h"
#include "curve.h"
#include "separatrix.h"
#include "gridzone.h"
#include "polsegm.h"

typedef struct EdgeSegment {
  DLListNode* tail;
  DLListNode* head;
  int n_size; //number of norm distribution 
  double* norm_dist;
  Curve* gridpoint_curve;
  struct EdgeSegment* next;
} EdgeSegment;

typedef struct {
  double xpt_r;
  double xpt_z;
  // the sequence of sep seg, the nomination rule.
  // 0 means which has interseciton with inner target.
  int index[4];
  EdgeSegment* edges[4]; 
} SepDistStr;

EdgeSegment* create_EdgeSegment();
SepDistStr* create_SepDistStr_from_sep(SeparatrixStr* sep);
void free_EdgeSegment(EdgeSegment* seg);
void free_SepDistStr(SepDistStr* sepdist);

void update_sn_SepDistStr_from_GridZoneInfo(SepDistStr* sepdist, GridZoneInfo* gzinfo);
void update_sn_SepDistStr_from_PolSegmsInfo(SepDistStr* sepdist, PolSegmsInfo* polseginfo);

#endif