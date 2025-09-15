#include "separatrix.h"

#include "ode.h"
#include <assert.h>
#include "mathbase.h"
#include <stdbool.h>

SeparatrixStr* init_separatrix_default(void)
{
  SeparatrixStr* sep = malloc(sizeof(SeparatrixStr));
  if (!sep)
  {
     fprintf(stderr, "Error: malloc failed in init_separatrix_default.\n");
    return NULL;
  }
  // initialize all value to zero
  sep->xpt_r=0.00;
  sep->xpt_z=0.00;
  sep->xpt_psi=0.00;
  sep->order = -1; // default value -1 means not use
  // the index is corespodng to the four line, later to sort the sequence.
  for (int i = 0; i<4; i++)
  {
    sep->index[i]=-1;
    sep->line_list[i] = NULL;
  }
  return sep;
}

void free_separatrix_default(SeparatrixStr* sep)
{
  if (!sep) return;
  for (int i = 0; i < 4; ++i) 
  {
    if (sep->line_list[i])
    {
      free_DLList(sep->line_list[i]);
    }
  }
  free(sep);
}
//generate the sepatrxi by a tracing method;
//interp is a 1D interpolator for calculate the psi in the X-point range
// func is a func that calculate the right side for a the ode solover, e.g. dy/dx = F, F is func.
// solver is a ode solver that based on different method, e.g. BRK45, euler etc.
// interp, func and solver should be constructed before use.
void generate_separatrix_bytracing(
  SeparatrixStr* sep, XPointInfo xpt, Equilibrium* equ, MagFieldTorSys* magfield,
  Interp1DFunction* interp, //bound to an interp interface
  ode_function* func,   // bound to a ode function
  ode_solver* solver  // dound to a ode solver
)
{
  int adjust=0; //use to expand the box which contain the X-point
  sep->xpt_psi=xpt->level;
  sep->xpt_r=xpt->centerX;
  sep->xpt_z=xpt->centerY;
  int cx1=xpt->cx1-adjust;
  int cx2=xpt->cx2+adjust;
  int cy1=xpt->cy1-adjust;
  int cy2=xpt->cy2+adjust;
  const double xpt_level = xpt->level;
  
  double level_tmp1;
  double level_tmp2;

  int x = cx1;
  int y = cy1;
  int ox, oy;

//TODO: check there are four points with psi value in the line of xpt rectangular.
  int nside=2*((cx2-cx1) + (cy2-cy1));
  

  int np=0;
  double second_p[4][2];
  int dir; //indicate the start point from X or Y side


//build the x_tmp,fx_tmp, dfdx_tmp will be used to update interpolator data.
  double x_tmp[2]={0};
  double fx_tmp[2]={0};
  double dfdx_tmp[2]={0}; 

//build the tracer;
  for(int i=0; i<nside;i++)
  {
    ox=x;oy=y;
    level_tmp1=equ->psi[ox][oy];
    if (y == cy1) {  // bottom edge
      if (x == cx2) {  // bottom-right corner → move up
          x_tmp[0] = equ->z[y];
          fx_tmp[0] = equ->psi[x][y];
          y++;
          x_tmp[1] = equ->z[y];
          fx_tmp[1] = equ->psi[x][y];
          dir=1;
          second_p[np][0] = equ->r[x];
          // printf("second_p[np][0]: %.15f\n",second_p[np][0]);
          // printf("second_p[np][1]: %.15f\n",second_p[np][1]);

      } else {  // move right
          x_tmp[0] = equ->r[x];
          fx_tmp[0] = equ->psi[x][y];
          x++;
          x_tmp[1] = equ->r[x];
          fx_tmp[1] = equ->psi[x][y];
          dir=0;
          second_p[np][1] = equ->z[y];
      }
    }
    else if (x == cx2) {  // right edge
      if (y == cy2) {  // top-right corner → move left
          x_tmp[0] = equ->r[x];
          fx_tmp[0] = equ->psi[x][y];
          x--;
          x_tmp[1] = equ->r[x];
          fx_tmp[1] = equ->psi[x][y];
          dir=0;
          second_p[np][1] = equ->z[y];
      } else {  // move up
          x_tmp[0] = equ->z[y];
          fx_tmp[0] = equ->psi[x][y];
          y++;
          x_tmp[1] = equ->z[y];
          fx_tmp[1] = equ->psi[x][y];
          dir=1;
          second_p[np][0] = equ->r[x];
      }
    }
    else if (y == cy2) {  // top edge
      if (x == cx1) {  // top-left corner → move down
          x_tmp[0] = equ->z[y];
          fx_tmp[0] = equ->psi[x][y];
          y--;
          x_tmp[1] = equ->z[y];
          fx_tmp[1] = equ->psi[x][y];
          dir=1;
          second_p[np][0] = equ->r[x];
      } else {  // move left
          x_tmp[0] = equ->r[x];
          fx_tmp[0] = equ->psi[x][y];
          x--;
          x_tmp[1] = equ->r[x];
          fx_tmp[1] = equ->psi[x][y];
          dir=0;
          second_p[np][1] = equ->z[y];
      }
    }
    else if (x == cx1) {  // left edge
      if (y == cy1) {  // bottom-left corner → move right
          x_tmp[0] = equ->r[x];
          fx_tmp[0] = equ->psi[x][y];
          x++;
          x_tmp[1] = equ->r[x];
          fx_tmp[1] = equ->psi[x][y];
          dir=0;
          second_p[np][1] = equ->z[y];

      } else {  // move down
          x_tmp[0] = equ->z[y];
          fx_tmp[0] = equ->psi[x][y];
          y--;
          x_tmp[1] = equ->z[y];
          fx_tmp[1] = equ->psi[x][y];
          dir=1;
          second_p[np][0] = equ->r[x];
      }
    }
    else {
      assert(0 && "Point is not on the boundary");
    }
    
    level_tmp2=equ->psi[x][y];

    //in some situation, the start point can be directly the point.
    double f1 = level_tmp1 - xpt_level;
    double f2 = level_tmp2 - xpt_level;
    double eps = EPSILON_12;
    
    bool crosses = (f1 * f2 <= 0.0) && (fabs(f1) > eps || fabs(f2) > eps);
  
    if (crosses)
    {
      //ensure x_tmp[0] < x_tmp[1] for interpolation
      if(x_tmp[0]>x_tmp[1])
      {
        swap_double(&x_tmp[0],&x_tmp[1]);
        swap_double(&fx_tmp[0],&fx_tmp[1]);
      }
      modify_cubicherm1D_data(interp->data, x_tmp, fx_tmp, NULL,2);
      Secant_Method(&second_p[np][dir],&xpt_level, interp);

      // printf("Second point: %.15f %.15f\n", second_p[np][0], second_p[np][1]);
      printf("Find the psi point in Xpt Rec:\n");
      printf("ox: %d x: %d oy: %d y: %d\n", ox, x, oy, y);
      printf("%.15f %.15f\n", second_p[np][0], second_p[np][1]);
      np++;
    }
    if(np==4) break; //find all four point
  }

  if(np!=4)
  {
    fprintf(stderr, "Unexpected error: don't find all the four points for separatrix lines.\n");
    exit(EXIT_FAILURE);
  }
  //Begin to trace;
  //check the direction, because the direction should far away from X-point

  double *next_p=malloc(2*sizeof(double));
  double *start_p=malloc(2*sizeof(double));
  char filename[100];
  for(int i=0;i<4;i++)
  {
    sprintf(filename, "sep_baseline%d", i);
    FILE *fp = fopen(filename, "w");
    start_p[0]=second_p[i][0];
    start_p[1]=second_p[i][1];
    double t=0.00;
    double step_size=0.1; //artifical number.

    //one step to check whether reverse the direction
    solver->next_step(step_size, &t, start_p, next_p, solver->solver_data, func);

    //check direction if the next point is in xpt rec.
    if((next_p[0]>equ->r[cx1]&&next_p[0]<equ->r[cx2]&&next_p[1]>equ->z[cy1]&&next_p[1]<equ->z[cy2])) 
    {
      printf("Change the direction for %d sep line.\n", i);
      for(int i=0;i<func->ndim;i++)
      {
        func->rescale[i]=-1.0*func->rescale[i];
      }
    }

    //print the four point for the Xpt rec for debug.
    // fprintf(fp,"%.15f %.15f\n",equ->r[cx1], equ->z[cy1]);
    // fprintf(fp,"%.15f %.15f\n",equ->r[cx1], equ->z[cy2]);
    // fprintf(fp,"%.15f %.15f\n",equ->r[cx2], equ->z[cy1]);
    // fprintf(fp,"%.15f %.15f\n",equ->r[cx2], equ->z[cy2]);

    //print the Xpoint and 2nd point for the sep line.
    fprintf(fp,"%.15f %.15f\n",xpt->centerX, xpt->centerY);
    fprintf(fp,"%.15f %.15f\n",start_p[0], start_p[1]);

    //Add X-point to separatix line
    sep->line_list[i]=create_DLListNode(xpt->centerX,xpt->centerY);
    //Add 2nd point to separatix line
    DLListNode* endnode=get_DLList_tailnode(sep->line_list[i]);
    add_DLListnode_at_tail(&endnode, start_p[0], start_p[1]);

    t=t+step_size;
    int boundary=10; //not reach the real boundary but [10:nx-10][10:ny-10]

    int tracing_counter=0;
    while(tracing_counter<MAX_NUM_TRACING)
    {
      start_p[0]=next_p[0];
      start_p[1]=next_p[1];
      t=t+step_size;
      solver->next_step(step_size, &t, start_p, next_p, solver->solver_data, func);
      tracing_counter+=1;
      if (next_p[0] < equ->r[boundary] || next_p[0] > equ->r[equ->nw - boundary] ||
          next_p[1] < equ->z[boundary] || next_p[1] > equ->z[equ->nh - boundary])
      {
        printf("Arrive the boundary!\n");
        break;
      }
      
      if((next_p[0] > equ->r[cx1] && next_p[0]<equ->r[cx2] &&
          next_p[1] > equ->z[cy1] && next_p[1]<equ->z[cy2]))
      {
        printf("Back to X-Point Rectangular!\n");
        //insert Xpoint to the sep lines which are the LCFS.
        add_DLListnode_at_tail(&endnode, xpt->centerX, xpt->centerY);
        fprintf(fp, "%.15f %.15f\n", xpt->centerX, xpt->centerY);
        break;
      }
      add_DLListnode_at_tail(&endnode, next_p[0], next_p[1]);
      fprintf(fp, "%.15f %.15f\n", next_p[0], next_p[1]);
    }
    if(tracing_counter==MAX_NUM_TRACING)
    {
      printf("Arrive the Maxium Tracing numbr: %d.\n", tracing_counter);
      printf("WARNING: Please DOUBLE CHECK the file %s.\n", filename);
    }


    //add the X-point again to make sure the sep lines for LCFS is closed.

    fclose(fp);
  }
  free(start_p);
  free(next_p); 
}