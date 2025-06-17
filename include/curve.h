#ifndef CURVE_H
#define CURVE_H

#include<stdlib.h>
#include<datastructure.h>

// Represents a single 2D curve with n_point points
typedef struct {
    size_t n_point;         // Number of (x, y) points in the curve
    double** points;      // Array of points [n_point][2], each point is [x, y]
} OldCurve;

// Represents a set of curves

// Allocate memory for a single curve with 'n_point' (x,y) pairs
OldCurve* create_oldcurve(size_t n_point);

// Free memory of a single curve
void free_oldcurve(OldCurve* c);

int copy_oldcurve(OldCurve* c1, OldCurve* c2);


typedef struct {
    double x, y;  // Coordinates of the point
    // double psi; // Optional: physical fields like flux, temperature, etc.
} CurvePoint;

/**
 * Represents a 2D curve with dynamic memory for storing points.
 */
typedef struct {
    size_t n_point;      // Number of points currently stored
    size_t capacity;     // Total allocated capacity
    CurvePoint* points;  // Dynamically allocated array of points
} Curve;

// Create a new curve with initial capacity (default to 128 if zero)
Curve* create_curve(size_t init_capacity);

// Free memory associated with the curve
void free_curve(Curve* c);

// Append a point (x, y) to the curve, automatically expanding memory
int add_last_point_curve(Curve* c, double x, double y);

// Set coordinates of the i-th point (must be within bounds)
int set_point_curve(Curve* c, size_t i, double x, double y);

int delete_last_point(Curve* c);

void write_curve(const char *filename, const Curve *c);

void print_curve(const Curve *c);


/*************************************************
 * Use for TwoDim Grid Generation
**************************************************/

//total length of the curve
double total_length_curve(const Curve *c);

//calcule the length until the nth points of a curve.
//Coresponding to 'long_CARRE'
double length_curve(const Curve *c, size_t n);

//Coresponding to 'indcrb_CARRE'
int indcrb_curve(const Curve* curve, const CurvePoint* point, double d);

//Coresponding to 'ruban_CARRE'
double ruban_curve(const Curve* curve, const CurvePoint* point, double d);



//calcaute the coordinates along the curve with the distance d;
void coordnates_in_curve(const Curve* curve, double d, CurvePoint* point);



#endif