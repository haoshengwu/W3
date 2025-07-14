#ifndef POLSEGM_H
#define POLSEGM_H
#include <stdio.h>
#include <stdlib.h>


typedef struct 
{
  char* name;
  int n_points;
  double* norm_dist;
}PolSegStr;


typedef struct 
{
  char* topo;
  int n_polsegms;
  PolSegStr** polsegments;
  // similar to GridZoneInfo, the following indicates
  // the polosegm corespond to which sep line part.
  int* xptidx; //which x-point
  int* seplineidx; //which lince
  int* segmidx; //the segment of the line, for multiple X-points

  //whethere reverse the direction to have the same with sep line 
  //Should always be 0 (false) for polsegms 
  int* reverse_segm;
}PolSegmsInfo;


PolSegmsInfo* create_PolSegmsInfo(int n_polsegms);
PolSegStr* create_PolSegStr(int n_points, char* name);

void free_PolSegStr(PolSegStr* polseg);
void free_PolSegmsInfo(PolSegmsInfo* polseginfo);

void write_PolSegStr(PolSegStr* polseg, FILE* fp);
void write_PolSegmsInfo(PolSegmsInfo* polseginfo, const char* filename);

void print_PolSegStr(PolSegStr* polseg);
void print_PolSegmsInfo(PolSegmsInfo* polseginfo);

PolSegmsInfo* read_PolSegmsInfo_from_file(const char* filename);
#endif