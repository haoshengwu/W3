#ifndef OPOINT_H
#define OPOINT_H


typedef struct {
  int cx1,cy1,cx2,cy2;
  double centerX,centerY;
  double psi; 
} OPointStr;

OPointStr* create_opoint();
void free_opoint(OPointStr* opoint);
#endif