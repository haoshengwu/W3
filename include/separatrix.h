#ifndef SEPARATRIX_H
#define SEPARATRIX_H

#include "datastructure.h"
#include "xpoint.h"
#include "equilibrium.h"
#include "magneticfield.h"
#include "ode.h"
#include "config.h"



typedef struct SeparatrixStr SeparatrixStr;
typedef struct SeparatrixOpt SeparatrixOpt;

// A Structure for Separatrix
typedef struct SeparatrixStr
{
    double xpt_r, xpt_z;             // X-point coordinates
    double xpt_psi;              // Magnetic flux at X-point
    // Indices defining separatrix topology, use to nominate the four sep lines
    //index[0] is the the sep line have intersection with Inner Target.
    //Then index[1],[2],[3] are begin in the counter-clockwise order.
    int index[4];            
    int order;               // 1st, 2nd, 3rd, or 4th X-point
    DLListNode* line_list[4]; // Pointers to 4 separatrix line segments
} SeparatrixStr;

// general interface for initialize, free, and generator separatrix.
typedef SeparatrixStr* (*Init_Separatrix_Fun)(void);
typedef void (*Free_Separatrix_Fun)(SeparatrixStr* sep);
typedef void (*Generate_Separatrix_Fun)(
    SeparatrixStr* sep, 
    _XPointInfo* xpt, 
    Equilibrium* equ,
    MagFieldTorSys* magfield,
    Interp1DFunction* interp, //bound to an interp interface
    void* func,   // bound to a ode function
    void* solver  // dound to a ode solver
);

typedef struct SeparatrixOpt
{
    Init_Separatrix_Fun init;       // Function to allocate/init separatrix
    Free_Separatrix_Fun free;       // Function to clean up
    Generate_Separatrix_Fun generate_sep;// Strategy function to generate separatrix
} SeparatrixOpt;

SeparatrixStr* init_separatrix_default(void);
void free_separatrix_default(SeparatrixStr* sep);
void generate_separatrix_bytracing(
  SeparatrixStr* sep, XPointInfo xpt, Equilibrium* equ, MagFieldTorSys* magfield,
  Interp1DFunction* interp, //bound to an interp interface
  ode_function* func,   // bound to a ode function
  ode_solver* solver // dound to a ode solver
);


#endif