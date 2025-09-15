#ifndef SEPDISTRIBUTION_H
#define SEPDISTRIBUTION_H

#include "datastructure.h"
#include "curve.h"
#include "separatrix.h"
#include "gridzoneinfo.h"
#include "polsegminfo.h"
#include "config.h"


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
  // 0 means the first 'spe line' which has interseciton with inner target.
  int index[4];
  EdgeSegment* edges[4]; 
  int order; //coresponding to order in sep. indicates the Xth of xpt
} SepDistStr;

EdgeSegment* create_EdgeSegment();
SepDistStr* create_SepDistStr_from_sep(SeparatrixStr* sep);
void free_EdgeSegment(EdgeSegment* seg);
void free_SepDistStr(SepDistStr* sepdist);

//using solgirdinfo to update sn sepdist
//Cut the sep lines
void update_sn_SepDistStr_from_GridZoneInfo(SepDistStr* sepdist, GridZoneInfo* gzinfo);

//using solgirdinfo to update sn sepdist
//update the normal distribution
void update_sn_SepDistStr_from_PolSegmsInfo(SepDistStr* sepdist, PolSegmsInfo* polseginfo);

// calcualte the gridpoint_curve for the SepDistStr itself
void update_SepDistStr_gridpoint_curve(SepDistStr* sepdist);

// create a gridpoint_curve according to head and norm_distribution, 
// n_norm_dist is the size of normal_distribution.

Curve* create_gridpoint_curve(DLListNode* head, double* norm_dist, int n_norm_dist);




#endif