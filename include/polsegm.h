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
}PolSegmsInfo;


PolSegmsInfo* create_PolSegmsInfo(int n_polsegms);
PolSegStr* create_PolSegStr(int n_points, char* name);

void free_PolSegStr(PolSegStr* polseg);
void free_PolSegmsInfo(PolSegmsInfo* polseginfo);

void write_PolSegStr(PolSegStr* polseg, FILE* fp);
void write_PolSegmsInfo(PolSegmsInfo* polseginfo, const char* filename);

PolSegmsInfo* read_PolSegmsInfo_from_file(const char* filename);
#endif