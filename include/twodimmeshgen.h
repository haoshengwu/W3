#ifndef TWODIMMESHGEN_H
#define TWODIMMESHGEN_H

#include "carrefunction.h"
#include "structuredgrid.h"
#include "ode.h"

void generate_CARRE_mesh(GridZone* GridZone, 
                         Equilibrium* equ,
                         ode_function* func, 
                         ode_solver* solver);






#endif