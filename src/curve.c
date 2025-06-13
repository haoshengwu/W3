#include "curve.h"
#include <stdio.h>

Curve* create_curve(size_t n_point)
{   
  if(!n_point)
  {
    fprintf(stderr, "Error: n_point is 0!\n");
    return NULL;
  }
  Curve* c = malloc(sizeof(Curve));
  if (!c)
  {
    fprintf(stderr, "Error: Memory allocation failed for a curve!\n");
    return NULL;
  } 
  c->n_point = n_point;
  c->points = allocate_2d_array(n_point,2);
  return c;
}

void free_curve(Curve* c) {
    if (!c) return;
    free_2d_array(c->points);  // Free row pointers
    free(c);         // Free the Curve struct
}

CurveSet* create_curveset(size_t n_curve, size_t n_point) {
    CurveSet* cs = malloc(sizeof(CurveSet));
    if (!cs) {
        fprintf(stderr, "Error: Failed to allocate CurveSet\n");
        return NULL;
    }

    cs->n_curve = n_curve;
    cs->curves = malloc(sizeof(Curve*) * n_curve);

    for (size_t i = 0; i < n_curve; i++) {
        cs->curves[i] = create_curve(n_point);  // Return Curve*
        if (!cs->curves[i]) {
            fprintf(stderr, "Error: Failed to create curve %zu\n", i);
            // Free previously allocated curves before exiting
            for (size_t j = 0; j < i; j++)
                free_curve(cs->curves[j]);
            free(cs->curves);
            free(cs);
            return NULL;
        }
    }

    return cs;
}

void free_curveset(CurveSet *cs)
{
  if (!cs) return;
  for(size_t i=0; i<cs->n_curve;i++)
  {
    free_curve(cs->curves[i]);
  }
  free(cs->curves);
  free(cs);
}

 int copy_curve(Curve* c1, Curve* c2)
{
    if (!c1 || !c2) {
        fprintf(stderr, "Null pointer passed to copy_curve.\n");
        return 1;
    }

    if (c1->n_point != c2->n_point) {
        fprintf(stderr, "The size of two curves is not the same!\n");
        return 1;
    }

    for (size_t i = 0; i < c1->n_point; i++) {
        c1->points[i][0] = c2->points[i][0];
        c1->points[i][1] = c2->points[i][1];
    }
    return 0;
}

