#ifndef SEPDISTRIBUTION_H
#define SEPDISTRIBUTION_H

#include "datastructure.h"
#include "curve.h"


typedef struct EdgeSegment {
  DLListNode* tail;
  DLListNode* head;
  int n_size;
  double* norm_dist;
  Curve* curve;
  struct EdgeSegment* next;
} EdgeSegment;

typedef struct XPoint {
  double r;
  double z;
  EdgeSegment* edges[4]; 
} SepDistStr;


#endif