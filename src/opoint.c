#include "opoint.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Allocate and initialize an OPointStr with default values:
// - Integers set to -1 (undefined index)
// - Doubles set to NaN (undefined physical values)
OPointStr* create_opoint()
{
    OPointStr* pt = (OPointStr*) malloc(sizeof(OPointStr));
    if (!pt) {
        fprintf(stderr, "Memory allocation failed in create_opoint_nan.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize integer fields to -1 (invalid/undefined index)
    pt->cx1 = -1;
    pt->cy1 = -1;
    pt->cx2 = -1;
    pt->cy2 = -1;

    // Initialize double fields to NaN (undefined physical values)
    pt->centerX = NAN;
    pt->centerY = NAN;
    pt->psi     = NAN;

    return pt;
}

// Free an OPointStr instance
void free_opoint(OPointStr* pt)
{
    if (pt) free(pt);
}