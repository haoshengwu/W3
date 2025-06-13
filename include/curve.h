#ifndef CURVE_H
#define CURVE_H

#include<stdlib.h>
#include<datastructure.h>

// Represents a single 2D curve with n_point points
typedef struct {
    size_t n_point;         // Number of (x, y) points in the curve
    double** points;      // Array of points [n_point][2], each point is [x, y]
} Curve;

// Represents a set of curves
typedef struct {
    size_t n_curve;      // Number of curves in the set
    Curve** curves;       // Array of Curve structures
} CurveSet;

// Allocate memory for a single curve with 'n_point' (x,y) pairs
Curve* create_curve(size_t n_point);

// Free memory of a single curve
void free_curve(Curve* c);

// Allocate memory for a curve set with 'n_curve' curves,
// each containing 'n_point' points
CurveSet* create_curveset(size_t n_curve, size_t n_point);

// Free memory of the entire curve set (all curves and their points)
void free_curveset(CurveSet* cs);

int copy_curve(Curve* c1, Curve* c2);

#endif