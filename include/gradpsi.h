#ifndef GRADPSI_H
#define GRADPSI_H

#include "datastructure.h"
#include "xpoint.h"
#include "equilibrium.h"
#include "magneticfield.h"
#include "ode.h"
#include "separatrix.h"
#include "opoint.h"

typedef struct {
    int nr;         // Number of grid points in the r direction
    int nz;         // Number of grid points in the z direction
    double *r;      // r-coordinate array (size nr)
    double *z;      // z-coordinate array (size nz)
    double **psi;   //psi(r,z)

    // First-order gradient of ψ (gradψ = [∂ψ/∂r, ∂ψ/∂z])
    double **gradpsi_r;      // ∂ψ/∂r
    double **gradpsi_z;      // ∂ψ/∂z

    // Second derivatives of gradψ with respect to r
    double **dgradpsi_r_dr;  // ∂(∂ψ/∂r)/∂r
    double **dgradpsi_z_dr;  // ∂(∂ψ/∂z)/∂r

    // Second derivatives of gradψ with respect to z
    double **dgradpsi_r_dz;  // ∂(∂ψ/∂r)/∂z
    double **dgradpsi_z_dz;  // ∂(∂ψ/∂z)/∂z

    // Mixed second derivatives
    double **d2gradpsi_r_drdz; // ∂²(∂ψ/∂r)/∂r∂z
    double **d2gradpsi_z_drdz; // ∂²(∂ψ/∂z)/∂r∂z

} GradPsiStr;

GradPsiStr* init_grad_psi(void);
void free_grad_psi(GradPsiStr *gradpsi);
void calc_grad_psi(Equilibrium *equ, GradPsiStr *gradpsi, diff2d_eval_fun func);
void write_grad_psi(GradPsiStr *gradpsi, const char* filename);




// A Structure for grad of psi lines
typedef struct 
{
    double xpt_r, xpt_z;             // X-point coordinates
    double xpt_psi;              // Magnetic flux at X-point
    int index[4];            // Indices defining separatrix topology
    int order;               // 1st, 2nd, 3rd, or 4th gradient psi lines
    DLListNode* line_list[4]; // Pointers to 4 grad psi line segments
} GradPsiLineStr;

GradPsiLineStr* init_gradpsiline_default(void);
void free_gradpsiline_default(GradPsiLineStr* gradpsi_lines);

void generate_gradpsiline_bytracing(
  GradPsiLineStr* gradpsi_lines,
  GradPsiStr* gradpsi,
  OPointStr* opoint,
  SeparatrixStr* sep,
  Interp1DFunction* interp, //NOT USE NOW, Can be update
  ode_function* func,   // bound to a ode function
  ode_solver* solver // dound to a ode solver
);

#endif