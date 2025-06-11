#include "gradpsi.h"
#include "mathbase.h"
#include "datastructure.h"
#include <math.h>


GradPsiStr* init_grad_psi(void)
{
    GradPsiStr* gradpsi=malloc(sizeof(GradPsiStr));
    gradpsi->nr=-1;
    gradpsi->nz=-1;
    gradpsi->r=NULL;
    gradpsi->z=NULL;
    gradpsi->psi=NULL;
    gradpsi->gradpsi_r=NULL;
    gradpsi->gradpsi_z=NULL;
    gradpsi->dgradpsi_r_dr=NULL;
    gradpsi->dgradpsi_r_dz=NULL;
    gradpsi->dgradpsi_z_dr=NULL;
    gradpsi->dgradpsi_z_dz=NULL;
    gradpsi->d2gradpsi_r_drdz=NULL;
    gradpsi->d2gradpsi_z_drdz=NULL;
    return gradpsi;
}

void free_grad_psi(GradPsiStr *gradpsi)
{
    if(!gradpsi) return;
    if(gradpsi->r) free(gradpsi->r);
    if(gradpsi->z) free(gradpsi->z);
    if(gradpsi->psi) free_2d_array(gradpsi->psi);
    if(gradpsi->gradpsi_r) free_2d_array(gradpsi->gradpsi_r);
    if(gradpsi->gradpsi_z) free_2d_array(gradpsi->gradpsi_z);
    if(gradpsi->dgradpsi_r_dr) free_2d_array(gradpsi->dgradpsi_r_dr);
    if(gradpsi->dgradpsi_r_dz) free_2d_array(gradpsi->dgradpsi_r_dz);
    if(gradpsi->dgradpsi_z_dr) free_2d_array(gradpsi->dgradpsi_z_dr);
    if(gradpsi->dgradpsi_z_dz) free_2d_array(gradpsi->dgradpsi_z_dz);
    if(gradpsi->d2gradpsi_r_drdz) free_2d_array(gradpsi->d2gradpsi_r_drdz);
    if(gradpsi->d2gradpsi_z_drdz) free_2d_array(gradpsi->d2gradpsi_z_drdz);
    free(gradpsi);
}

void calc_grad_psi(Equilibrium *equ, GradPsiStr *gradpsi, diff2d_eval_fun func)
{
// func(nx, r[], ny, z[], f[nx][ny], dfdx[nx][ny], dfdy[nx][ny])
  int nr = equ->nw;
  int nz = equ->nh;
  gradpsi->nr = nr;
  gradpsi->nz = nz;
  gradpsi->psi=allocate_2d_array(nr,nz);
  gradpsi->r = malloc(gradpsi->nr*sizeof(double));
  gradpsi->z = malloc(gradpsi->nz*sizeof(double));

  for(int i=0; i<nr; i++)
  {
    gradpsi->r[i]=equ->r[i];
  }

  for(int j=0; j<nz; j++)
  {
    gradpsi->z[j]=equ->z[j];
  }

  for(int i=0; i<nr; i++)
  {
    for (int j=0; j<nz; j++)
    {
      gradpsi->psi[i][j]=equ->psi[i][j];
    }
  }

  gradpsi->gradpsi_r=allocate_2d_array(nr,nz);
  gradpsi->gradpsi_z=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_r_dr=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_r_dz=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_z_dr=allocate_2d_array(nr,nz);
  gradpsi->dgradpsi_z_dz=allocate_2d_array(nr,nz);
  gradpsi->d2gradpsi_r_drdz=allocate_2d_array(nr,nz);
  gradpsi->d2gradpsi_z_drdz=allocate_2d_array(nr,nz);
  func(nr, gradpsi->r, nz, gradpsi->z, gradpsi->psi, gradpsi->gradpsi_r, gradpsi->gradpsi_z);
  func(nr, gradpsi->r, nz, gradpsi->z,gradpsi->gradpsi_r, gradpsi->dgradpsi_r_dr, gradpsi->dgradpsi_r_dz);
  func(nr, gradpsi->r, nz, gradpsi->z,gradpsi->gradpsi_z, gradpsi->dgradpsi_z_dr, gradpsi->dgradpsi_z_dz);

  for (int i = 0; i < nr; i++) 
  {
    for (int j = 0; j < nz; j++) 
    {
      int ixp = min(nr-1, i+1);
      int ixm = max(0,     i-1);
      int iyp = min(nz-1, j+1);
      int iym = max(0,     j-1);

      double dr = gradpsi->r[ixp] - gradpsi->r[ixm];
      double dz = gradpsi->z[iyp] - gradpsi->z[iym];
      double scale = 1.0 / (dr * dz);

      gradpsi->d2gradpsi_r_drdz[i][j] = (
          gradpsi->gradpsi_r[ixp][iyp] - gradpsi->gradpsi_r[ixm][iyp]
        - gradpsi->gradpsi_r[ixp][iym] + gradpsi->gradpsi_r[ixm][iym]
      ) * scale;

      gradpsi->d2gradpsi_z_drdz[i][j] = (
        gradpsi->gradpsi_z[ixp][iyp] - gradpsi->gradpsi_z[ixm][iyp]
        - gradpsi->gradpsi_z[ixp][iym] + gradpsi->gradpsi_z[ixm][iym]
      ) * scale;
    }
  }
  printf("Finish calculating gradpsi.\n");
}

void write_grad_psi(GradPsiStr *gradpsi, const char* filename)
{
  if(!gradpsi||!filename)
  {
    fprintf(stderr,"Empty input parameters for write_grad_psi.\n");
    exit(EXIT_FAILURE);
  }
  FILE* fp;
  fp=fopen(filename, "w");
  int nr=gradpsi->nr;
  int nz=gradpsi->nz;
  fprintf(fp, "# r z gradpsi_x gradpsi_y\n");
  for(int i=0; i<nr; i++)
  {
    for(int j=0; j<nz; j++)
    {
      fprintf(fp, "%lf %lf %lf %lf\n",
              gradpsi->r[i],gradpsi->z[j],
              gradpsi->gradpsi_r[i][j], gradpsi->gradpsi_z[i][j]);
    }
  }
  fclose(fp);
  printf("Successfully writing gradpsi to %s.\n",filename);
}

GradPsiLineStr* init_gradpsiline_default(void)
{
  GradPsiLineStr* gradpsi_lines = malloc(sizeof(GradPsiLineStr));
  if (!gradpsi_lines)
  {
     fprintf(stderr, "Error: malloc failed in init_separatrix_default.\n");
    return NULL;
  }
  //initial values
  gradpsi_lines->xpt_r=0.0;
  gradpsi_lines->xpt_z=0.0;
  gradpsi_lines->order=-1;
  for(int i=0; i<4; i++)
  {
    gradpsi_lines->index[i]=-1;
    gradpsi_lines->line_list[i]=NULL;
  }
  return gradpsi_lines;
}

void free_gradpsiline_default(GradPsiLineStr* gradpsi_lines)
{
  if(!gradpsi_lines) return;
  for (int i = 0; i < 4; ++i) 
  {
    if (gradpsi_lines->line_list[i])
    {
      free_DLList(gradpsi_lines->line_list[i]);
    }
  }
  free(gradpsi_lines);
}

//first line (x2-x1,y2-y1) second line(x3-x2, y3-y2)
// 0 means the same direction, 1 means not.
static int same_direction(double x1, double y1, double x2, double y2, double x3, double y3)
{
  double line1[2];
  double line2[2];
  line1[0]=x2-x1; line1[1]=y2-y1;
  line2[0]=x3-x2; line2[1]=y3-y2;
  double direction = line1[0]*line2[0]+line1[1]*line2[1];

  if(direction>0)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}
static void normalize(double v[2]) {
    double len = sqrt(v[0]*v[0] + v[1]*v[1]);
    if (len != 0.0) {
        v[0] /= len;
        v[1] /= len;
    }
}

void compute_angle_bisector_point(double x0, double y0,
                                  double x1, double y1,
                                  double x2, double y2,
                                  double* x3, double* y3)
{
    // Vectors v1 = P1 - P0, v2 = P2 - P0
    double v1[2] = { x1 - x0, y1 - y0 };
    double v2[2] = { x2 - x0, y2 - y0 };

    // Compute lengths
    double len1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
    double len2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
    double length_avg = 0.5 * (len1 + len2);

    // Normalize vectors
    normalize(v1);
    normalize(v2);

    // Bisector direction = normalized sum of unit vectors
    double bisector[2] = { v1[0] + v2[0], v1[1] + v2[1] };
    normalize(bisector);

    // Final point = origin + average length * bisector direction
    *x3 = x0 + length_avg * bisector[0];
    *y3 = y0 + length_avg * bisector[1];
}

void generate_gradpsiline_bytracing(
  GradPsiLineStr* gradpsi_lines,
  GradPsiStr* gradpsi,
  OPointStr* opoint,
  SeparatrixStr* sep,
  Interp1DFunction* interp, //NOT USE NOW
  ode_function* func,   // bound to a ode function
  ode_solver* solver // dound to a ode solver
)
{
  //prepare the X-point as the first point of the four lines
  double xpt_r = sep->line_list[0]->r;
  double xpt_z = sep->line_list[0]->z;
  gradpsi_lines->xpt_psi = sep->xpt_psi;

  printf("DEBUG xpt_r xpt_z\n");
  printf("%lf %lf\n", xpt_r, xpt_z);

  //use for tracing
  double *next_p=malloc(2*sizeof(double));
  double *start_p=malloc(2*sizeof(double));
  //file to record
  char filename[100];

  printf("Begin to generate gardpsi lines by tracing.\n");
  for(int i=0; i<4; i++)
  {
    sprintf(filename, "gradpsiline_baseline%d", i);
    FILE *fp = fopen(filename, "w");
    double t=0.00;
    double step_size=0.01; //artifical number.

    //the second point is decided by the 2nd points which contain the gradpsi line
    //                 \    :     /
    //                  \   :    /
    //            ...... X-point......gradpsi line
    //                  /   :   \
    //                 /    :    \
    
    DLListNode* sec_node1=sep->line_list[i]->next;
    DLListNode* sec_node2=sep->line_list[(i+1)%4]->next;

    //current solution is not solid, gradpsi line is average of the adjunct two sep lines.
    double length=1;
    compute_angle_bisector_point(xpt_r, xpt_z,
                             sec_node1->r, sec_node1->z,
                             sec_node2->r, sec_node2->z,
                             &start_p[0], &start_p[1]);

    //one step to check whether reverse the direction
    solver->next_step(step_size, &t, start_p, next_p, solver->solver_data, func);
    if(!same_direction(xpt_r, xpt_z, start_p[0], start_p[1], next_p[0], next_p[1])==0)
    {
      printf("DEBUG change direction for gradpsi line%d \n", i);
      for(int i=0;i<func->ndim;i++)
      {
        func->rescale[i]=-1.0*func->rescale[i];
      }
    }

    gradpsi_lines->line_list[i]=create_DLListNode(xpt_r, xpt_z);
    DLListNode* endnode=get_DLList_endnode(sep->line_list[i]);
    insert_DLList_at_end(&endnode, start_p[0], start_p[1]);

    t=t+step_size;
    int boundary=10; //psi range for stopping tracing.
    double opoint_delta=0.1; //opoint range for stopping tracing.

    while(1)
    {
      start_p[0]=next_p[0];
      start_p[1]=next_p[1];
      t=t+step_size;
      solver->next_step(step_size, &t, start_p, next_p, solver->solver_data, func);

      if (next_p[0] < gradpsi->r[boundary] || next_p[0] > gradpsi->r[gradpsi->nr - boundary] ||
          next_p[1] < gradpsi->z[boundary] || next_p[1] > gradpsi->z[gradpsi->nz - boundary])
      {
        printf("Arrive the boundary!\n");
        break;
      }

      if((next_p[0] > opoint->centerX-opoint_delta &&
          next_p[0] < opoint->centerX+opoint_delta &&
          next_p[1] > opoint->centerY-opoint_delta &&
          next_p[1] < opoint->centerY+opoint_delta))
      {
        printf("Back to X-Point Rectangular!\n");
        break;
      }

      //Forced termination because gradpsi may be too high
      double dis = sqrt(pow(next_p[0]-start_p[0],2)+pow(next_p[1]-start_p[1],2));
      if(dis < 1.0E-6)
      {
        printf("Gradient of psi is too small, stop tracing.\n");
        break;
      }

      insert_DLList_at_end(&endnode, next_p[0], next_p[1]);
      fprintf(fp, "%.15f %.15f\n", next_p[0], next_p[1]);
    }
    printf("Finish tracing of gradpsi line %d \n", i);
    fclose(fp);
  }

  free(start_p);
  free(next_p); 
}
