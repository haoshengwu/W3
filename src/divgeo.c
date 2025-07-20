#include "divgeo.h"
#include "datastructure.h"
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "gridzoneinfo.h"
#include "utils.h"
#include "polsegminfo.h"

#define MAX_ITER_DG 1000
#define EPSILON_DG 5.0E-7

#define TGUARD1 0.1
#define TGUARD2 0.1
#define PASMIN 1.0E-3



DGRadRegion* create_DGRadRegion(int n_level, const char* name)
{
  DGRadRegion* radregion = malloc(sizeof(DGRadRegion));
  if (!radregion) {
    fprintf(stderr, "Error: failed to allocate memory for DGRadRegion\n");
    return NULL;
  }

  // Allocate memory for name
  radregion->name = malloc(strlen(name) + 1);  // +1 for null terminator
  if (!radregion->name) 
  {
    fprintf(stderr, "Error: failed to allocate memory for radregion name\n");
    free(radregion);
    return NULL;
  }
  strcpy(radregion->name, name);  // Now safe


  radregion->n_level = n_level;
  radregion->level = malloc(sizeof(double) * n_level);
  if (!radregion->level) {
    fprintf(stderr, "Error: failed to allocate memory for radregion levels\n");
    free(radregion);
    return NULL;
  }

  return radregion;
}

void free_DGRadRegion(DGRadRegion* radregion)
{
  if (!radregion) return;

  if (radregion->name) 
  {
    free(radregion->name);
    radregion->name = NULL;
  }

  if (radregion->level) 
  {
    free(radregion->level);
    radregion->level = NULL;
  }

  free(radregion);
}

DGPolZone* create_DGPolZone(int n_points, const char* name)
{
  if(!name)
  {
    fprintf(stderr, "Error: The name is NULL\n");
    return NULL;
  }
  DGPolZone* polzone = malloc(sizeof(DGPolZone));
  if (!polzone) 
  {
    fprintf(stderr, "Error: failed to allocate memory for DGPolZone\n");
    return NULL;
  }

  polzone->name = malloc(strlen(name) + 1);
  if (!polzone->name) 
  {
    fprintf(stderr, "Error: failed to allocate memory for polzone name\n");
    free(polzone);
    return NULL;
  }
  strcpy(polzone->name, name);

  polzone->n_points = n_points;
  polzone->norm_dist = malloc(sizeof(double) * n_points);
  if (!polzone->n_points) 
  {
    fprintf(stderr, "Error: failed to allocate memory for polzone levels\n");
    free(polzone->name);
    free(polzone);
    return NULL;
  }

  return polzone;
}

void free_DGPolZone(DGPolZone* polzone)
{
  if (!polzone) return;

  if (polzone->name) {
    free(polzone->name);
    polzone->name = NULL;
  }

  if (polzone->norm_dist) {
    free(polzone->norm_dist);
    polzone->norm_dist = NULL;
  }
  free(polzone);
}

DivGeoTrg* create_dgtrg(void)
{
  DivGeoTrg* trg=malloc(sizeof(DivGeoTrg));
  if (!trg) 
  {
    fprintf(stderr, "Error: failed to allocate memory for DivGeoTrg\n");
    return NULL;
  }
  trg->topo==NULL;
  trg->target_curves=NULL;
  trg->regions=NULL;
  trg->zones=NULL;

  trg->dltr1=NULL;
  trg->dltrn=NULL;
  trg->npr=NULL;

  trg->dltp1=NULL;
  trg->dltpn=NULL;
  trg->nptseg=NULL;

  return trg;
}

//not use now
static void strip_newline(char* line) {
    line[strcspn(line, "\r\n")] = '\0';
}

//not use now
static int is_keyword(const char* line, const char* keyword) {
    return strncmp(line, keyword, strlen(keyword)) == 0;
}

static DGRadRegion* read_region(FILE* fp, int n_level, const char* name)
{
    char line[256];
    DGRadRegion* radregion = create_DGRadRegion(n_level, name);
    radregion->level[0]=NAN;
    fgets(line, sizeof(line), fp);//skip the line which is 'level'
    for(int i=1; i<n_level; i++)
    {
      fgets(line, sizeof(line), fp);
      line[strcspn(line, "\r\n")] = '\0';
      sscanf(line, "%lf",&(radregion->level[i]));
    }
    printf("%s\n",radregion->name);
    printf("n_level: %d\n",n_level);
    for(int i=0; i<n_level; i++)
    {
      printf("%lf\n", radregion->level[i]);
    }
    return radregion;
}

static DGPolZone* read_zone(FILE* fp, int n_points, const char* name)
{
    char line[256];

    DGPolZone* polzone = create_DGPolZone(n_points, name);
    polzone->norm_dist[0]=0.00;
    polzone->norm_dist[n_points-1]=1.00;

    fgets(line, sizeof(line), fp);//skip the line which is 'points'

    for(int i=1; i<n_points-1; i++)
    {
      fgets(line, sizeof(line), fp);
      line[strcspn(line, "\r\n")] = '\0';
      sscanf(line, "%lf", &(polzone->norm_dist[i]));
    }
    printf("%s", polzone->name);
    printf("n_points: %d\n",n_points);

    for(int i=0; i<n_points; i++)
    {
      printf("%lf\n", polzone->norm_dist[i]);
    }
    return polzone;
}

static OldCurve* read_curve(FILE* fp, int* n_curve_point)
{
    char line[256];
    int n_point = 0;
    long pos;
    while (fgets(line, sizeof(line), fp)) 
    {
      line[strcspn(line, "\r\n")] = '\0';
      // printf("DEBUG read_curve line: %s\n", line);
      if (strncmp(line, "line", 4)==0 ) pos = ftell(fp);
      else if (strncmp(line, "target", 6)==0 || strncmp(line, "region", 6)==0) break;
      else n_point++;
    }
    //printf("Number of points for target cure: %d\n",n_point);
    *n_curve_point=n_point;
    fseek(fp, pos, SEEK_SET);

    OldCurve* curve = create_oldcurve(n_point);
    for (int i = 0; i < n_point; i++) 
    {
      fgets(line, sizeof(line), fp);
      sscanf(line, "%lf , %lf", &(curve->points[i][0]), &(curve->points[i][1]));
      //convert from mm to m!!!!
      curve->points[i][0]=curve->points[i][0]/1000.0;
      curve->points[i][1]=curve->points[i][1]/1000.0;
      printf("%lf %lf\n",curve->points[i][0],curve->points[i][1]);
    }
    return curve;
}

static double* read_doubles(FILE* fp, int n) {
    if (n <= 0) return NULL;

    double* values = malloc(n * sizeof(double));
    if (!values) {
        fprintf(stderr, "Failed to allocate memory for doubles.\n");
        return NULL;
    }

    char line[256];
    int count = 0;
    for(int i=0;i<n;i++)
    {
      fgets(line, sizeof(line), fp);
      line[strcspn(line, "\r\n")] = '\0';
      sscanf(line,"%lf",&(values[i]));
      printf("%lf\n",values[i]);
    }

    return values;
}

static int* read_ints(FILE* fp, int n) {
    int* values = malloc(n * sizeof(int));
    if (!values) {
        fprintf(stderr, "Failed to allocate memory for int array.\n");
        return NULL;
    }

    char line[256];
    int count = 0;
    while (count < n && fgets(line, sizeof(line), fp)) {
        line[strcspn(line, "\r\n")] = '\0';
        values[count++] = atoi(line);
    }

    if (count < n) {
        fprintf(stderr, "Expected %d integers, but got %d\n", n, count);
        free(values);
        return NULL;
    }

    return values;
}

int load_dgtrg_from_file(DivGeoTrg* trg, const char* filename)
{
    if (!trg || !filename) {
        fprintf(stderr, "Input parameters are invalid!\n");
        return 1;
    }

    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        return 1;
    }

    // First pass: count and collect data for preallocation
    char line[256];
    int  count_zones = 0, count_regions = 0, count_targets = 0;
    int* temp_npr = NULL;
    int* temp_nptseg = NULL;

    while (fgets(line, sizeof(line), fp)) {
        line[strcspn(line, "\r\n")] = '\0';
        if (strncmp(line, "zone", 4) == 0) count_zones++;
        else if (strncmp(line, "region", 6) == 0) count_regions++;
        else if (strncmp(line, "target", 6) == 0) count_targets++;
        else if (strcmp(line, "npr") == 0) {
            if (temp_npr) free(temp_npr);
            temp_npr = read_ints(fp, count_regions);
        }
        else if (strcmp(line, "nptseg") == 0) {
            if (temp_nptseg) free(temp_nptseg);
            temp_nptseg = read_ints(fp, count_zones);
        }
    }
    printf("Target curve number: %d\n", count_targets);
    printf("Region number: %d\n", count_regions);
    printf("Zone number: %d\n", count_zones);

    for(int i=0; i<count_regions; i++)
    {
      printf("npr: region %d %d\n", i, temp_npr[i]);
    }

    for(int i=0; i<count_zones; i++)
    {
      printf("nptseg: zone %d %d\n", i, temp_nptseg[i]);
    }

    rewind(fp);

    int cur_target = 0, cur_region = 0, cur_zone = 0;

    trg->n_target = count_targets;
    trg->n_region = count_regions;
    trg->n_zone = count_zones;

    trg->n_target_curve = count_targets > 0 ? malloc(count_targets * sizeof(int)) : NULL;
    trg->target_curves = count_targets > 0 ? malloc(count_targets * sizeof(OldCurve*)) : NULL;
    trg->regions = count_regions > 0 ? malloc(count_regions * sizeof(DGRadRegion*)) : NULL;
    trg->zones = count_zones > 0 ? malloc(count_zones * sizeof(DGPolZone*)) : NULL;

    trg->dltr1 = trg->dltrn = NULL;
    trg->npr = NULL;
    trg->dltp1 = trg->dltpn = NULL;
    trg->nptseg = NULL;
    trg->topo = NULL;

    // Use already read npr and nptseg arrays
    trg->npr = temp_npr;
    trg->nptseg = temp_nptseg;
    char temp_name[32];
    while (fgets(line, sizeof(line), fp)) 
    {
      line[strcspn(line, "\r\n")] = '\0';
//      printf("DEBUG line: %s\n",line);
//      printf("DEBUG strncmp: %d\n",strncmp(line, "# topo", 6));
      
      if (strncmp(line, "# topo", 6) == 0) 
      {
//        printf("DEBUG in the loop.\n");
        char* topo_name = line + 6;
        while (*topo_name == ' ' || *topo_name == '\t') topo_name++;
        topo_name[strcspn(topo_name, "\r\n")] = '\0';
        trg->topo = malloc(strlen(topo_name) + 1);
        if (!trg->topo) { fprintf(stderr, "Memory allocation failed for topo name.\n"); goto fail;}
        strcpy(trg->topo, topo_name);
        printf("Topology: %s\n",trg->topo);
        continue;
      }
      
      if (line[0] == '#' || strlen(line) == 0) continue;

      //Read target curves
      if (strncmp(line, "target", 6) == 0) 
      {
        printf("Read target curve %d\n", cur_target);
        printf("%s\n", line);
        trg->target_curves[cur_target] = read_curve(fp,&trg->n_target_curve[cur_target]);
        printf("Number of points in the target curve: %d\n", trg->n_target_curve[cur_target]);
        printf("Finish reading target curve %d\n", cur_target);
        cur_target++;
        continue;
      }

      //Read radial regions
      if (strncmp(line, "region", 6) == 0)
      {
        snprintf(temp_name, sizeof(temp_name), "region%d", cur_region);
        // printf("%s\n",temp_name);
        trg->regions[cur_region] = read_region(fp, trg->npr[cur_region],temp_name);
        printf("Finish reading region %d\n", cur_region);
        cur_region++;
        continue;
      }

      //Read poloidal zones
      if (strncmp(line, "zone", 4) == 0) 
      {
        snprintf(temp_name, sizeof(temp_name), "zone%d\n", cur_zone);
        trg->zones[cur_zone] = read_zone(fp, trg->nptseg[cur_zone],temp_name);
        printf("Finish reading zone %d\n", cur_zone);
        cur_zone++;
        continue;
      }
      if (strcmp(line, "dltr1") == 0) 
      { 
        printf("%s:\n",line);
        trg->dltr1 = read_doubles(fp, trg->n_zone); 
        continue; 
      }
      if (strcmp(line, "dltrn") == 0) 
      { 
        printf("%s:\n",line);
        trg->dltrn = read_doubles(fp, trg->n_zone);
        continue; 
      }
      if (strcmp(line, "pntrat") == 0) 
      {
        fgets(line, sizeof(line), fp);
        trg->pntrat = atof(line);
        continue;
      }
      if (strcmp(line, "dltp1") == 0) 
      { 
        printf("%s:\n",line);
        trg->dltp1 = read_doubles(fp, trg->n_region); 
        continue; 
      }
      if (strcmp(line, "dltpn") == 0) 
      { 
        printf("%s:\n",line);
        trg->dltpn = read_doubles(fp, trg->n_region); 
        continue; 
      }
      if (strncmp(line, "clstruct", 8) == 0) break;
    }
    fclose(fp);
    return 0;

fail:
    if (fp) fclose(fp);
    free_dgtrg(trg);
    return 1;
}

TargetDLListCurve* create_target_curve_from_dgtrg(DivGeoTrg* trg, int n)
{
  TargetDLListCurve* tgt_cur=create_target_curve();
  if (!tgt_cur) 
  {
    fprintf(stderr, "Failed to allocate memmory for target_curve");
    exit(EXIT_FAILURE);
  }
  char name[32];
  snprintf(name, sizeof(name), "TargetCurve_%d", n);
  change_name_target_curve(tgt_cur, name);
  int n_point=trg->n_target_curve[n];
  //  printf("Target OldCurve Name: %s\n", tgt_cur->name);
  for(int i=0; i<n_point;i++)
  {
    // printf("DEBUG r=%.4f z=%.4f\n", r, z);
    double r = trg->target_curves[n]->points[i][0];
    double z = trg->target_curves[n]->points[i][1];
    // printf("DEBUG i=%d\n", i);
    // printf("DEBUG r=%.4f z=%.4f\n", r, z);
    add_point_target_curve(tgt_cur, r, z);
  }
  return tgt_cur;
}

void free_dgtrg(DivGeoTrg* trg)
{
  if (!trg) return;

  if (trg->topo) 
  {
    free(trg->topo);
    trg->topo = NULL;
  }

  if (trg->target_curves) 
  {
    for (int i = 0; i < trg->n_target; i++) 
    {
      if (trg->target_curves[i]) 
      {
        free_oldcurve(trg->target_curves[i]);
        trg->target_curves[i] = NULL;
      }
    }
    free(trg->target_curves);
    trg->target_curves = NULL;
  }

  if (trg->regions) 
  {
    for (int i = 0; i < trg->n_region; i++) 
    {
      if (trg->regions[i]) 
      {
        free_DGRadRegion(trg->regions[i]);
        trg->regions[i] = NULL;
      }
    }
    free(trg->regions);
    trg->regions = NULL;
  }

  if (trg->zones) 
  {
    for (int i = 0; i < trg->n_zone; i++) 
    {
      if (trg->zones[i]) 
      {
        free_DGPolZone(trg->zones[i]);
        trg->zones[i] = NULL;
      }
    }
    free(trg->zones);
    trg->zones = NULL;
  }

  free(trg->n_target_curve); trg->n_target_curve=NULL;
  free(trg->dltr1); trg->dltr1 = NULL;
  free(trg->dltrn); trg->dltrn = NULL;
  free(trg->npr);   trg->npr = NULL;

  free(trg->dltp1);  trg->dltp1 = NULL;
  free(trg->dltpn);  trg->dltpn = NULL;
  free(trg->nptseg); trg->nptseg = NULL;

  free(trg);
}


// along the curve which start from DDLListNode* head, 
//calcuate the corespoding r, z for the psi value based on the equilibirum and interp2d1f function;
// the total number of psi is n_psi.
//We assume the psi is monotic along the curve line.
//interp2d1f is used to calculate the psi at any position.
//We also assume the target cureve is linear.
//DLListNode* head is TargetDLListCurve head
//int skipfirst=1 means skip the first psi, which is correspoding to sep.
//int skiplast=1 means skip the last psi, which is correspoding to 2nd sep for Snwoflake/multiple X-points.

static void cal_points_from_psi(double *psi,double *r, double *z, int n,
                             DLListNode* head,
                             Equilibrium* equ,
                             interpl_2D_1f interp2d1f,
                             int skipfirst, int skiplast)
{

  double psi_head, psi_tail;
  DLListNode* end = get_DLList_tailnode(head);
  interp2d1f(head->r,head->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_head, NULL, NULL, NULL);
  interp2d1f(end->r,end->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_tail, NULL, NULL, NULL);
  printf("DEBUG psi_head %.15f psi_tail %.15f\n", psi_head, psi_tail);

  if((psi[0]>psi_head-1.0E-6&&psi[0]<psi_tail+1.0E-6) 
      ||(psi[0]<psi_head+1.0E-6&&psi[0]>psi_tail-1.0E-6))
  {
    printf("psi[0] %.15fis in the range psi_head %.15f psi_tail %.15f with delta %.15f\n", psi[0], psi_head, psi_tail ,EPSILON_DG);
  }
  else
  {
    fprintf(stderr, "The psi[0]> %.15f is out of range psi_head %.15f psi_tail %.15f with delta %.15f\n!", psi[0], psi_head, psi_tail,EPSILON_DG);
    exit(EXIT_FAILURE);
  }

  if((psi[n-1]>psi_head-1.0E-6&&psi[n-1]<psi_tail+1.0E-6) 
      ||(psi[n-1]<psi_head+1.0E-6&&psi[n-1]>psi_tail-1.0E-6))
  {
    printf("psi[n-1] %.15f is in the range psi_head %.15f psi_tail %.15f with delta %.15f\n", psi[n-1], psi_head, psi_tail,EPSILON_DG);
  }
  else
  {
    fprintf(stderr, "psi[n-1] %.15f is out of range psi_head %.15f psi_tail %.15f with delta %.15f\n!", psi[n-1], psi_head, psi_tail,EPSILON_DG);
    exit(EXIT_FAILURE);
  }
  if(fabs(psi[1]-psi[0])<1.0E-5)
  {
    fprintf(stderr,"first magnetic surface is too close to separatrix, psi[1]-psi[0] < 1.0E-5. \n");
    exit(EXIT_FAILURE);
  }

  for(int i=0; i<n; i++)
  {
    if((i==0&&skipfirst)||(i==n-1&&skiplast))
    {
      continue;
    }
    DLListNode* current = head;
    printf("DEBUG create start point: %d\n", i);
    double psi_prev;
    double psi_curr;
    double psi_next;

    while(current->next!=NULL)
    {
      interp2d1f(current->r,current->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_prev, NULL, NULL, NULL);
      interp2d1f(current->next->r,current->next->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_curr, NULL, NULL, NULL);
      if((psi[i]>psi_prev-EPSILON_DG&&psi[i]<psi_curr+EPSILON_DG)
         ||(psi[i]<psi_prev+EPSILON_DG&&psi[i]>psi_curr-EPSILON_DG)) 
      {
        break; //Find the correct point
      }
      else current=current->next;
    }
    //use Secant_Method to calcute the position and assume linear between heat and heat->next

    double t_prev=0.0;
    double t_curr=1.0;
    double t_next;
    double psi_target=psi[i];
    for (int iter = 0; iter < MAX_ITER_DG; ++iter) 
    {
      double denom = psi_curr - psi_prev;
      if (fabs(denom) < EPSILON_DG) 
      {
        printf("Warning: function difference too small in Secant iteration.\n");
        break;
      }
      t_next = t_curr - (psi_curr - psi_target) * (t_curr - t_prev) / denom;
      double r_next = current->r + (current->next->r-current->r)*t_next;
      double z_next = current->z + (current->next->z-current->z)*t_next;
      interp2d1f(r_next,z_next, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_next, NULL, NULL, NULL);
      if (fabs(psi_next - psi_target) < EPSILON_DG) 
      {
        r[i] = r_next;
        z[i] = z_next;
        break;
      }
      t_prev=t_curr;
      psi_prev = psi_curr;
      t_curr=t_next;
      psi_curr=psi_next;
    }
  }
}


//0 means monotonic psi else 1.
static int check_curve_monotonic_psi(DLListNode* head,
                                     Equilibrium* equ,
                                     interpl_2D_1f interp2d1f)
{
  bool increasing = true;
  bool decreasing = true;
  DLListNode* newhead=head;
  while(newhead->next!=NULL)
  {
    double psi;
    interp2d1f(newhead->r,newhead->z,equ->nw,equ->r, equ->nh, equ->z, equ->psi, &psi, NULL, NULL,NULL);
    double nextpsi;
    interp2d1f(newhead->next->r,newhead->next->z,equ->nw,equ->r, equ->nh, equ->z, equ->psi, &nextpsi, NULL, NULL,NULL);
    if (fabs(psi-nextpsi)<EPSILON_DG)
    {
      fprintf(stderr,"The psi values are too close!\n");
      exit(EXIT_FAILURE);
    }
    else if(psi>nextpsi) 
    {
      increasing=false;
    }
    else if (psi<nextpsi)
    {
      decreasing=false;
    }

    newhead = newhead->next;
  }
  if(increasing||decreasing)
  {
    return 0; // monotonic
  }
  else
  {
    return 1; // not monotonic
  }
}

static double* compute_differences(const double* arr, int n)
{
    if (!arr || n < 2)
    {
        fprintf(stderr, "Invalid input for compute_differences.\n");
        exit(EXIT_FAILURE);
    }

    double* diff = (double*)malloc((n - 1) * sizeof(double));
    if (!diff)
    {
        fprintf(stderr, "Memory allocation failed in compute_differences.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n - 1; ++i)
    {
        diff[i] = arr[i + 1] - arr[i];
    }

    return diff;
}

static void reverse_array(double* arr, size_t n)
{
    if (!arr || n < 2) return;

    size_t i = 0, j = n - 1;
    while (i < j)
    {
        double tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
        i++;
        j--;
    }
}

static void scale_array(double* arr, size_t n, double factor)
{
    if (!arr || n == 0) return;

    for (size_t i = 0; i < n; ++i)
    {
        arr[i] *= factor;
    }
}

static double* connect_arrays(const double* arr1, size_t n1,
                              const double* arr2, size_t n2,
                              const double* arr3, size_t n3,
                              const double* arr4, size_t n4)
{
    size_t total = n1 + n2 + n3 + n4;

    if ((n1 && !arr1) || (n2 && !arr2) || (n3 && !arr3) || (n4 && !arr4))
    {
        fprintf(stderr, "Null array passed to connect_arrays.\n");
        exit(EXIT_FAILURE);
    }

    double* result = malloc(total * sizeof(double));
    if (!result)
    {
        fprintf(stderr, "Memory allocation failed in connect_arrays.\n");
        exit(EXIT_FAILURE);
    }

    size_t offset = 0;
    if (n1) { for (size_t i = 0; i < n1; ++i) result[offset++] = arr1[i]; }
    if (n2) { for (size_t i = 0; i < n2; ++i) result[offset++] = arr2[i]; }
    if (n3) { for (size_t i = 0; i < n3; ++i) result[offset++] = arr3[i]; }
    if (n4) { for (size_t i = 0; i < n4; ++i) result[offset++] = arr4[i]; }

    return result;
}

static double* recover_array_from_diff(const double* diff, int n, double first_value)
{
    if (!diff || n < 2)
    {
        fprintf(stderr, "Invalid input for recover_array_from_diff.\n");
        exit(EXIT_FAILURE);
    }

    double* arr = malloc((n+1) * sizeof(double));
    if (!arr)
    {
        fprintf(stderr, "Memory allocation failed in recover_array_from_diff.\n");
        exit(EXIT_FAILURE);
    }

    arr[0] = first_value;
    for (int i = 1; i < n+1; ++i)
    {
        arr[i] = arr[i - 1] + diff[i - 1];
    }

    return arr;
}

static void change_name(char **name, const char *str)
{
    if (*name) {
        free(*name);
        *name = NULL;
    }

    *name = malloc(strlen(str) + 1);
    if (*name != NULL) {
        strcpy(*name, str);
    } else {
        fprintf(stderr, "Failed to allocate memory\n");
        exit(EXIT_FAILURE);
    }
}

//Not used yet
// static DLListNode* cal_pol_points(DLListNode* head,
//                                   const double*     pol_dis,
//                                   int               n)
// {
//     if (!head || !pol_dis || n <= 0)
//     {
//       fprintf(stderr,"Empty input for cal_pol_points.\n");
//       exit(EXIT_FAILURE);
//     }
//     // 1. Compute total arc length of the original curve
//     double total_len = total_length_DLList(head);
//    DLListNode* new_head = NULL;
//     DLListNode* new_tail = NULL;
//     // 2. Interpolate points at each normalized distance
//     for (int k = 0; k < n; ++k) {
//         const double target_d = pol_dis[k] * total_len;
//         const DLListNode* seg = head;
//         double accum = 0.0;
//         DLListNode* new_node = NULL;
//         // Walk through the original curve to find the segment where target_d lies
//         while (seg && seg->next) {
//             double dx = seg->next->r - seg->r;
//             double dy = seg->next->z - seg->z;
//             double seg_len = sqrt(dx*dx + dy*dy);
//             // Check if target_d is within this segment, with floating point tolerance
//             if (accum - 1.0E-10 <= target_d && target_d <= accum + seg_len +  1.0E-10)
//             {
//                 double remain = target_d - accum;
//                 if (fabs(remain) <= 1.0E-10) {
//                     // Point lies at the start of the segment
//                     new_node = create_DLListNode(seg->r, seg->z);
//                 } else if (fabs(remain - seg_len) <= 1.0E-10) {
//                     // Point lies at the end of the segment
//                     new_node = create_DLListNode(seg->next->r, seg->next->z);
//                 } else {
//                     // Linear interpolation within the segment
//                     double ratio = remain / seg_len;
//                     double r_int = seg->r + ratio * dx;
//                     double z_int = seg->z + ratio * dy;
//                     new_node = create_DLListNode(r_int, z_int);
//                 }
//                 break;
//             }
//             accum += seg_len;
//             seg = seg->next;
//         }
//         // If no segment matched due to floating point rounding, use the last point as fallback
//         if (!new_node && seg) {
//             const DLListNode* tail = seg->next ? seg->next : seg;
//             new_node = create_DLListNode(tail->r, tail->z);
//         }
//         // Append the new node to the output linked list
//         if (!new_head) {
//             new_head = new_tail = new_node;
//         } else {
//             new_tail->next = new_node;
//             new_node->prev = new_tail;
//             new_tail = new_node;
//         }
//     }
//     if(pol_dis[0]<1.0E-10)
//     {
//       new_head->r=head->r;
//       new_head->z=head->z;
//     }
//     if((pol_dis[n-1]-1.0)<1.0E-10)
//     {
//       DLListNode* tail = get_DLList_tailnode(head);
//       new_tail->r=tail->r;
//       new_tail->z=tail->z;
//     }
//     return new_head;
// }

void write_sn_gridzoneinfo_from_dgtrg(DivGeoTrg* trg, Equilibrium* equ, SeparatrixStr* sep, GradPsiLineStr* gradpsilines)
{
  
  printf("Begion to create inputfils from DivGeo trg file.\n");

  printf("DEBUG topology: %s\n", trg->topo);
  //check the trg is for SNL topology.
  if(strcmp(trg->topo,"SNL"))
  {
    fprintf(stderr, "This is not SNL topology.\n");
    exit(EXIT_FAILURE);
  }
  if(trg->n_region!=3 || trg->n_region!=3)
  {
    fprintf(stderr, "The number of Zone and Regions are not consistent with SNL.\n");
  }

/***********************************************
*    STEP0 create GridZoneInfo
***********************************************/
  GridZoneInfo* solgzinfo=allocate_GridZoneInfo();
  GridZoneInfo* pfrgzinfo=allocate_GridZoneInfo();
  GridZoneInfo* coregzinfo=allocate_GridZoneInfo();

  update_GridZoneInfo_from_dgtrg(solgzinfo, trg, 0);
  update_GridZoneInfo_from_dgtrg(pfrgzinfo, trg, 1);
  update_GridZoneInfo_from_dgtrg(coregzinfo, trg, 2);

/***********************************************
*    STEP1 create TargetDLListCurve
***********************************************/
  int idx_inner_tgt = 1; //start from 0, the second target is the inner target

  TargetDLListCurve* inner_tgt_curve=create_target_curve_from_dgtrg(trg, idx_inner_tgt);
  printf("Inner Target name: %s\n", inner_tgt_curve->name);

  int idx_outer_tgt = 0; //start from 0, the second target is the inner target

  TargetDLListCurve* outer_tgt_curve=create_target_curve_from_dgtrg(trg, idx_outer_tgt);
  
  
/***********************************************
*    STEP2 Sort index for sep and gradpsi lines
***********************************************/
  sort_sep_gradpsiline_by_targetcurve(inner_tgt_curve, sep, gradpsilines);

/*************************************************
*    STEP3 Create the target curve for SOL PFR and CORE
*************************************************/
  double itsct_r_inner, itsct_z_inner;
  
  //calcute the intersection point between sep and the inner target
  insert_intersections_DLList(sep->line_list[sep->index[0]],
                              inner_tgt_curve->head,
                              &itsct_r_inner, &itsct_z_inner);

  update_number_target_curve(inner_tgt_curve);
  write_DLList(inner_tgt_curve->head,"inner_targetcurve");

  update_GridZoneInfo_start_curve(solgzinfo,inner_tgt_curve);
  update_GridZoneInfo_start_curve(pfrgzinfo,inner_tgt_curve);

  //cut the separatrix which intersect with inner target
  cut_DLList_after_point(sep->line_list[sep->index[0]],itsct_r_inner, itsct_z_inner);

  double itsct_r_outer, itsct_z_outer;
  //calcute the intersection point between sep and the outer target
  insert_intersections_DLList(sep->line_list[sep->index[1]],
                              outer_tgt_curve->head,
                              &itsct_r_outer, &itsct_z_outer);
  update_number_target_curve(outer_tgt_curve);
  write_DLList(outer_tgt_curve->head,"outer_targetcurve");

/****************************************************************************
* Update ene curve for SOL and PFR
* CORE GridZoneInfo will be updated later by the start points itself
****************************************************************************/

  update_GridZoneInfo_end_curve(solgzinfo,outer_tgt_curve);
  update_GridZoneInfo_end_curve(pfrgzinfo,outer_tgt_curve);

  //cut the separatrix
  cut_DLList_after_point(sep->line_list[sep->index[1]],itsct_r_outer, itsct_z_outer);
  
  //create the curve which will used to store the SOL radregion
  TargetDLListCurve* sol_tgt_curve=create_target_curve();

  //create the SOL curve and the original inner target became for the PFR curve
  split_intersections_target_curve(inner_tgt_curve, itsct_r_inner, itsct_z_inner, sol_tgt_curve);
  
  //reverse the PFR curve which start form sep
  reverse_DLList_in_target_curve(inner_tgt_curve);

  //change the name
  change_name_target_curve(inner_tgt_curve, "PFR");
  change_name_target_curve(sol_tgt_curve, "SOL");

  //create the core curve based on gradpsilines. 2 is the index[2]. The third line!
  TargetDLListCurve* core_curve=create_core_curve_from_gradpsilines(gradpsilines, 2);
  change_name_target_curve(core_curve, "CORE");

  printf("%s target curve has %d points\n", inner_tgt_curve->name, inner_tgt_curve->n);
  printf("%s target curve has %d points\n", sol_tgt_curve->name, sol_tgt_curve->n);
  printf("%s target curve has %d points\n", core_curve->name, core_curve->n);

  // write_DLList(sep->line_list[sep->index[0]],"debug_sep1");
  // write_DLList(inner_tgt_curve->head,"debug_pfr");
  // write_DLList(sol_tgt_curve->head,"debug_sol")
  
  //update end curve
  

/*************************************************
*    STEP4 Calculate the start points used for trcing
*************************************************/
  if(check_curve_monotonic_psi(sol_tgt_curve->head, equ, cubicherm2d1f))
  {
    fprintf(stderr, "Target Curve in not monotonic psi.\n");
    exit(EXIT_FAILURE);
  }
  if(check_curve_monotonic_psi(inner_tgt_curve->head, equ, cubicherm2d1f))
  {
    fprintf(stderr, "Target Curve in not monotonic psi.\n");
    exit(EXIT_FAILURE);
  }

  double *r_tmp=malloc(trg->regions[0]->n_level*sizeof(double));
  double *z_tmp=malloc(trg->regions[0]->n_level*sizeof(double));

  r_tmp[0]=itsct_r_inner;
  z_tmp[0]=itsct_z_inner;

  cal_points_from_psi(trg->regions[0]->level, r_tmp, z_tmp, trg->regions[0]->n_level,
                      sol_tgt_curve->head, equ, cubicherm2d1f, 1,0);
  write_array2f(r_tmp,z_tmp,  trg->regions[0]->n_level, "SOL_start_point");
  update_GridZoneInfo_start_points(solgzinfo, r_tmp, z_tmp, solgzinfo->nr);
  write_array2f(solgzinfo->start_point_R, solgzinfo->start_point_Z, 
                       solgzinfo->nr,"SOLGR_RZ");


  cal_points_from_psi(trg->regions[1]->level, r_tmp, z_tmp, trg->regions[1]->n_level,
                      inner_tgt_curve->head, equ, cubicherm2d1f, 1,0);
  write_array2f(r_tmp,z_tmp,  trg->regions[0]->n_level, "PFR_start_point");
  update_GridZoneInfo_start_points(pfrgzinfo, r_tmp, z_tmp, pfrgzinfo->nr);
  write_array2f(pfrgzinfo->start_point_R, pfrgzinfo->start_point_Z, 
                       pfrgzinfo->nr,"PFRGR_RZ");


  r_tmp[0]=gradpsilines->line_list[0]->r;
  z_tmp[0]=gradpsilines->line_list[0]->z;
  cal_points_from_psi(trg->regions[2]->level, r_tmp, z_tmp, trg->regions[2]->n_level,
                      core_curve->head, equ, cubicherm2d1f, 1,0);
  write_array2f(r_tmp,z_tmp,  trg->regions[0]->n_level, "CORE_start_point");
  update_GridZoneInfo_start_points(coregzinfo, r_tmp, z_tmp, coregzinfo->nr);
  write_array2f(coregzinfo->start_point_R, coregzinfo->start_point_Z, 
                       coregzinfo->nr,"COREGR_RZ");

  //Use start points to update coregzinfo end curve;
  update_COREGridZoneInfo_start_curve(coregzinfo);
  update_COREGridZoneInfo_end_curve(coregzinfo);

/********************************************************************************
*  NO NEED anymore  
*  STEP5 Calculate the first poloidal points and normalize poloidal distribution
********************************************************************************/
// //The main idea here is read the lines and then connect them. 
// //And also for the distribution. because in DG, the distribution given in separated line

//   //1. FOR SOL radregion
//   //for DLList curve
//   DLListNode* SOL_bnd1=copy_DLList(sep->line_list[sep->index[0]]);
//   reverse_DLList(&SOL_bnd1);
//   // write_DLList(SOL_bnd1,"DEBUG_SOL_BND1");
//   //for distribution
//   int tmp_np1=trg->zones[0]->n_points;
//   double len_bnd1=total_length_DLList(SOL_bnd1);
//   double* dtbt_bnd1=compute_differences(trg->zones[0]->norm_dist, tmp_np1);
//   scale_array(dtbt_bnd1, tmp_np1-1, len_bnd1);
//   reverse_array(dtbt_bnd1, tmp_np1-1);

//   //for poloidal points
//   DLListNode* SOL_pol_p1=cal_pol_points(sep->line_list[sep->index[0]],trg->zones[0]->norm_dist, tmp_np1);
//   reverse_DLList(&SOL_pol_p1);
//   // write_DLList(SOL_pol_p1,"DEBUG_SOL_pol_p1");



//   //BECAREFUL, the TWO LCFS are not the exactly the same
//   //please use the same one in the follow

//   //for DLList curve
//   int nLCFS=3; //here we use the 3(The FOUTH becasue start from 0)
//   DLListNode* SOL_bnd2=copy_DLList(sep->line_list[sep->index[nLCFS]]);
//   reverse_DLList(&SOL_bnd2);

//   //for distribution
//   int tmp_np2=trg->zones[2]->n_points;
//   double len_bnd2=total_length_DLList(SOL_bnd2);
//   double* dtbt_bnd2=compute_differences(trg->zones[2]->norm_dist, tmp_np2);
//   scale_array(dtbt_bnd2, tmp_np2-1, len_bnd2);
//   reverse_array(dtbt_bnd2, tmp_np2-1);

//   //for poloidal points
//   DLListNode* SOL_pol_p2=cal_pol_points(sep->line_list[sep->index[nLCFS]],trg->zones[2]->norm_dist, tmp_np2);
//   reverse_DLList(&SOL_pol_p2);


//   //for DLList curve
//   DLListNode* SOL_bnd3=copy_DLList(sep->line_list[sep->index[1]]);
//   //for distribution
//   int tmp_np3=trg->zones[1]->n_points;
//   double* dtbt_bnd3=compute_differences(trg->zones[1]->norm_dist, tmp_np3);
//   double len_bnd3=total_length_DLList(SOL_bnd3);
//   scale_array(dtbt_bnd3, tmp_np3-1, len_bnd3);
//   //for poloidal points
//   DLListNode* SOL_pol_p3=cal_pol_points(sep->line_list[sep->index[1]],trg->zones[1]->norm_dist, tmp_np3);


//   //Final DLList curve
//   connect_DLList(SOL_bnd1, &SOL_bnd2,1);
//   connect_DLList(SOL_bnd1, &SOL_bnd3,1);
//   write_DLList(SOL_bnd1,"SOL_bnd");

//   //Final distribution
//   double* dtbt_SOL=connect_arrays(dtbt_bnd1,tmp_np1-1,dtbt_bnd2, tmp_np2-1,
//                                   dtbt_bnd3,tmp_np3-1, NULL, 0);
//   int npsol=tmp_np1-1+tmp_np2-1+tmp_np3-1;
//   scale_array(dtbt_SOL, npsol, 1/(len_bnd1+len_bnd2+len_bnd3));
//   double* dtbt_norm_SOL=recover_array_from_diff(dtbt_SOL, npsol, 0.0);
//   write_array(dtbt_norm_SOL, npsol+1, "SOL_norm_pol_dtbt");

//   //Final poloidal points
//   connect_DLList(SOL_pol_p1, &SOL_pol_p2,1);
//   connect_DLList(SOL_pol_p1, &SOL_pol_p3,1);
//   write_DLList(SOL_pol_p1,"SOL_pol_point");

//   //2. FOR PFR radregion

//   DLListNode* PFR_bnd1=copy_DLList(sep->line_list[sep->index[0]]);
//   reverse_DLList(&PFR_bnd1);
//   DLListNode* PFR_bnd2=copy_DLList(sep->line_list[sep->index[1]]);
//   connect_DLList(PFR_bnd1, &PFR_bnd2,1);
//   write_DLList(PFR_bnd1,"PFR_bnd");
  

//   //Final distribution
//   double* dtbt_PFR=connect_arrays(dtbt_bnd1,tmp_np1-1,dtbt_bnd3,tmp_np3-1,NULL,0, NULL, 0);
//   int nppfr=tmp_np1-1+tmp_np3-1;
//   scale_array(dtbt_PFR, nppfr, 1/(len_bnd1+len_bnd3));
//   double* dtbt_norm_PFR=recover_array_from_diff(dtbt_PFR, nppfr, 0.0);
//   write_array(dtbt_norm_PFR, nppfr+1, "PFR_norm_pol_dtbt");

//   //poloidal points
//   DLListNode* PFR_pol_p1=cal_pol_points(sep->line_list[sep->index[0]],trg->zones[0]->norm_dist, tmp_np1);
//   reverse_DLList(&PFR_pol_p1);
//   DLListNode* PFR_pol_p2=cal_pol_points(sep->line_list[sep->index[1]],trg->zones[1]->norm_dist, tmp_np1);
//   connect_DLList(PFR_pol_p1, &PFR_pol_p2,1);
//   write_DLList(PFR_pol_p1,"PFR_pol_point");



//   //3. FOR COREE radregion
//   DLListNode* CORE_bnd1=copy_DLList(sep->line_list[sep->index[nLCFS]]);
//   write_DLList(CORE_bnd1,"CORE_bnd");

//   //Final distribution
//   int npcore=tmp_np2-1;
//   double* dtbt_CORE=connect_arrays(dtbt_bnd2,tmp_np2-1,NULL,0,NULL,0, NULL, 0);
//   scale_array(dtbt_CORE, npcore, 1/(len_bnd2));
//   double* dtbt_norm_CORE=recover_array_from_diff(dtbt_CORE, npcore, 0.0);
//   write_array(dtbt_norm_CORE, npcore+1, "CORE_norm_pol_dtbt");

//   //Poloidal points
//   DLListNode* CORE_pol_p1=cal_pol_points(sep->line_list[sep->index[nLCFS]],trg->zones[2]->norm_dist, tmp_np2);
//   write_DLList(CORE_pol_p1,"CORE_pol_point");



//   update_GridZone_pol_norm_distrb(solgzinfo, dtbt_norm_SOL, npsol+1);
//   update_GridZone_pol_norm_distrb(pfrgzinfo, dtbt_norm_PFR, nppfr+1);
//   update_GridZone_pol_norm_distrb(coregzinfo, dtbt_norm_CORE, npcore+1);

//   update_GridZone_first_pol_points(solgzinfo, SOL_pol_p1);
//   update_GridZone_first_pol_points(pfrgzinfo, PFR_pol_p1);
//   update_GridZone_first_pol_points(coregzinfo, CORE_pol_p1);

/********************************************************************************
*    STEP6 Write the input based on GridZoneInfo
********************************************************************************/
  write_GridZoneInfo(solgzinfo,"gridzoneinfo_SOL");
  write_GridZoneInfo(pfrgzinfo,"gridzoneinfo_PFR");
  write_GridZoneInfo(coregzinfo,"gridzoneinfo_CORE");

/*************************************************
*    Free Memory 
*************************************************/
  free(r_tmp);
  free(z_tmp);
  free_target_curve(inner_tgt_curve);
  free_target_curve(outer_tgt_curve);

  free_target_curve(sol_tgt_curve);
  free_target_curve(core_curve);
  free_GridZoneInfo(&solgzinfo);
  free_GridZoneInfo(&pfrgzinfo);
  free_GridZoneInfo(&coregzinfo);

  // free_DLList(SOL_bnd1);
  // free_DLList(CORE_bnd1);
  // free_DLList(PFR_bnd1);

  // free_DLList(SOL_pol_p1);
  // free_DLList(PFR_pol_p1);
  // free_DLList(CORE_pol_p1);


  // free(dtbt_bnd1);
  // free(dtbt_bnd2);
  // free(dtbt_bnd3);
  // free(dtbt_SOL);
  // free(dtbt_norm_SOL);
  // free(dtbt_PFR);
  // free(dtbt_norm_PFR);
  // free(dtbt_CORE);
  // free(dtbt_norm_CORE);
}


void update_GridZoneInfo_from_dgtrg(GridZoneInfo* gridzoneinfo, DivGeoTrg* trg, int index)
{
  if (!gridzoneinfo || !trg) 
  {
    fprintf(stderr, "Empty input for update_GridZoneInfo_from_dgtrg.\n");
    exit(EXIT_FAILURE);
  }

  change_name(&(gridzoneinfo->topo), trg->topo);

  double guard_top;
  double guard_end;
 
  if (strcmp(trg->topo, "SNL") == 0) 
  {
    printf("Update GridZoneInfo in SNL topo.\n");

    if (index == 0) 
    {
      change_name(&(gridzoneinfo->name), "SOLGRIDZONEINFO");
      guard_top=TGUARD1;
      guard_end=TGUARD2;
    } 
    else if (index == 1) 
    {
      change_name(&(gridzoneinfo->name), "PFRGRIDZONEINFO");
      guard_top=TGUARD1;
      guard_end=TGUARD2;
    } 
    else if (index == 2) 
    {
      change_name(&(gridzoneinfo->name), "COREGRIDZONEINFO");
      guard_top=0.0;
      guard_end=0.0;
    } 
    else 
    {
      fprintf(stderr, "WRONG index in update_GridZoneInfo_from_dgtrg.\n");
      exit(EXIT_FAILURE);
    }
  }

  int nr = trg->npr[index];
  gridzoneinfo->nr = nr;

  gridzoneinfo->guard_start = malloc(nr * sizeof(double));
  gridzoneinfo->guard_end = malloc(nr * sizeof(double));
  gridzoneinfo->pasmin = malloc(nr * sizeof(double));

  if (!gridzoneinfo->guard_start || !gridzoneinfo->guard_end || !gridzoneinfo->pasmin) 
  {
    fprintf(stderr, "Failed to allocate memory for guard/pasmin.\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < nr; i++) 
  {
    gridzoneinfo->guard_start[i] = guard_top;
    gridzoneinfo->guard_end[i] = guard_end;
    gridzoneinfo->pasmin[i] = PASMIN;
  }

  // === Setup poloidal segments ===
  if (strcmp(trg->topo, "SNL") == 0) 
  {
    if (index == 0) 
    {  // SOL
      gridzoneinfo->n_polsegm1 = 3;
      gridzoneinfo->xptidx1 = malloc(3 * sizeof(int));
      gridzoneinfo->seplineidx1 = malloc(3 * sizeof(int));
      gridzoneinfo->segmidx1 = malloc(3 * sizeof(double*));
      gridzoneinfo->reverse_segm1 = malloc(3 * sizeof(int));
      gridzoneinfo->xptidx1[0]=0;
      gridzoneinfo->xptidx1[1]=0;
      gridzoneinfo->xptidx1[2]=0;
      gridzoneinfo->seplineidx1[0]=0; 
      gridzoneinfo->seplineidx1[1]=2;
      gridzoneinfo->seplineidx1[2]=1;
      gridzoneinfo->segmidx1[0]=0;
      gridzoneinfo->segmidx1[1]=0;
      gridzoneinfo->segmidx1[2]=0;
      gridzoneinfo->reverse_segm1[0]=1;
      gridzoneinfo->reverse_segm1[1]=1;
      gridzoneinfo->reverse_segm1[2]=0;

    }
    else if (index == 1) 
    {  // PFR
      gridzoneinfo->n_polsegm1 = 2;
      gridzoneinfo->xptidx1 = malloc(2 * sizeof(int));
      gridzoneinfo->seplineidx1 = malloc(2 * sizeof(int));
      gridzoneinfo->segmidx1 = malloc(2 * sizeof(double*));
      gridzoneinfo->reverse_segm1 = malloc(2 * sizeof(int));
      gridzoneinfo->xptidx1[0]=0;
      gridzoneinfo->xptidx1[1]=0;
      gridzoneinfo->seplineidx1[0]=0; 
      gridzoneinfo->seplineidx1[1]=1;
      gridzoneinfo->segmidx1[0]=0;
      gridzoneinfo->segmidx1[1]=0;
      gridzoneinfo->reverse_segm1[0]=1;
      gridzoneinfo->reverse_segm1[1]=0;

    }
    else if (index == 2) 
    {  // CORE
      gridzoneinfo->n_polsegm1 = 1;
      gridzoneinfo->xptidx1 = malloc( sizeof(int));
      gridzoneinfo->seplineidx1 = malloc(sizeof(int));
      gridzoneinfo->segmidx1 = malloc(sizeof(double*));
      gridzoneinfo->reverse_segm1 = malloc(sizeof(int));
      gridzoneinfo->xptidx1[0]=0;
      gridzoneinfo->seplineidx1[0]=2; 
      gridzoneinfo->segmidx1[0]=0;
      gridzoneinfo->reverse_segm1[0]=1;
    }
    else {
      fprintf(stderr, "WRONG index in update_GridZoneInfo_from_dgtrg.\n");
      exit(EXIT_FAILURE);
    }
  }
}

void update_GridZoneInfo_start_points(GridZoneInfo* gridzoneinfo, double* r, double* z, int n_point)
{
  if(!gridzoneinfo||!r||!z)
  {
    fprintf(stderr,"Empty input for update_GridZone_start_points.\n");
    exit(EXIT_FAILURE);
  }
  if(gridzoneinfo->nr!=n_point)
  {
    fprintf(stderr,"The points number is not consistent with GridZoneInfo\n");
    exit(EXIT_FAILURE);
  }
  //Check is there is already existing values
  if(gridzoneinfo->start_point_R || gridzoneinfo->start_point_Z)
  {
    fprintf(stderr,"There are existing values of R and Z.\n");
    exit(EXIT_FAILURE);
  }
  
  printf("DEBUGE n_point: %d\n", n_point);
  gridzoneinfo->start_point_R = malloc(n_point*sizeof(double));
  gridzoneinfo->start_point_Z = malloc(n_point*sizeof(double));

  if(!gridzoneinfo->start_point_R || !gridzoneinfo->start_point_Z)
  {
    fprintf(stderr,"Failed to alloc memmory for R and Z.\n");
    exit(EXIT_FAILURE);
  }

  for(int i=0; i<n_point; i++)
  {
    gridzoneinfo->start_point_R[i]=r[i];
    gridzoneinfo->start_point_Z[i]=z[i];
  }
}

void update_GridZoneInfo_start_curve(GridZoneInfo* gridzoneinfo, const TargetDLListCurve* tgt_cur)
{
    if (!gridzoneinfo || !tgt_cur || !tgt_cur->head || tgt_cur->n < 2)
    {
        fprintf(stderr, "Invalid input to update_GridZone_end_curve.\n");
        exit(EXIT_FAILURE);
    }

    if (gridzoneinfo->start_curve)
    {
        fprintf(stderr, "WARNING: There is alread an end_cure.\n");
        exit(EXIT_FAILURE);
    }

    int n = tgt_cur->n;
    //gridzoneinfo->end_curve = create_oldcurve(n);
    gridzoneinfo->start_curve = create_curve(n);
    if (!gridzoneinfo->start_curve)
    {
        fprintf(stderr, "Memory allocation failed for end_curve.\n");
        exit(EXIT_FAILURE);
    }

    DLListNode* node = tgt_cur->head;
    for (int i = 0; i < n; ++i)
    {
        if (!node)  // enough node
        {
            fprintf(stderr, "tgt_cur->n is larger than actual linked list length.\n");
            exit(EXIT_FAILURE);
        }
        add_last_point_curve(gridzoneinfo->start_curve, node->r,node->z);
        node = node->next;
    }

    if (node != NULL)
    {
        fprintf(stderr, "tgt_cur->n is smaller than actual linked list length.\n");
        exit(EXIT_FAILURE);
    }
}

void update_GridZoneInfo_end_curve(GridZoneInfo* gridzoneinfo, const TargetDLListCurve* tgt_cur)
{
    if (!gridzoneinfo || !tgt_cur || !tgt_cur->head || tgt_cur->n < 2)
    {
        fprintf(stderr, "Invalid input to update_GridZone_end_curve.\n");
        exit(EXIT_FAILURE);
    }

    if (gridzoneinfo->end_curve)
    {
        fprintf(stderr, "WARNING: There is alread an end_cure.\n");
        exit(EXIT_FAILURE);
    }

    int n = tgt_cur->n;
    //gridzoneinfo->end_curve = create_oldcurve(n);
    gridzoneinfo->end_curve = create_curve(n);
    if (!gridzoneinfo->end_curve)
    {
        fprintf(stderr, "Memory allocation failed for end_curve.\n");
        exit(EXIT_FAILURE);
    }

    DLListNode* node = tgt_cur->head;
    for (int i = 0; i < n; ++i)
    {
        if (!node)  // enough node
        {
            fprintf(stderr, "tgt_cur->n is larger than actual linked list length.\n");
            exit(EXIT_FAILURE);
        }
        add_last_point_curve(gridzoneinfo->end_curve, node->r,node->z);
        node = node->next;
    }

    if (node != NULL)
    {
        fprintf(stderr, "tgt_cur->n is smaller than actual linked list length.\n");
        exit(EXIT_FAILURE);
    }
}

void update_COREGridZoneInfo_end_curve(GridZoneInfo* gridzoneinfo)
{
    if (!gridzoneinfo || !gridzoneinfo->start_point_R 
        ||!gridzoneinfo->start_point_Z || gridzoneinfo->nr < 2) 
    {
        fprintf(stderr, "Empty or invalid input to update_COREGridZone_end_curve.\n");
        exit(EXIT_FAILURE);
    }
  if(gridzoneinfo->end_curve)
  {
    fprintf(stderr,"There is already end curve for COREgridzone\n");
  }
  
  int n = gridzoneinfo->nr;
  //gridzoneinfo->end_curve = create_oldcurve(n);  // Create first, then assign points
  gridzoneinfo->end_curve = create_curve(n);  // Create first, then assign points


  if (!gridzoneinfo->end_curve)
  {
    fprintf(stderr, "Memory allocation failed for end_curve.\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0; i<n; i++)
  {
    add_last_point_curve(gridzoneinfo->end_curve, 
                         gridzoneinfo->start_point_R[i], 
                         gridzoneinfo->start_point_Z[i]);
  }
}

void update_COREGridZoneInfo_start_curve(GridZoneInfo* gridzoneinfo)
{
    if (!gridzoneinfo || !gridzoneinfo->start_point_R 
        ||!gridzoneinfo->start_point_Z || gridzoneinfo->nr < 2) 
    {
        fprintf(stderr, "Empty or invalid input to update_COREGridZone_end_curve.\n");
        exit(EXIT_FAILURE);
    }
  if(gridzoneinfo->start_curve)
  {
    fprintf(stderr,"There is already end curve for COREgridzone\n");
  }
  
  int n = gridzoneinfo->nr;
  gridzoneinfo->start_curve = create_curve(n);  // Create first, then assign points

  if (!gridzoneinfo->start_curve)
  {
    fprintf(stderr, "Memory allocation failed for end_curve.\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0; i<n; i++)
  {
    add_last_point_curve(gridzoneinfo->start_curve, 
                         gridzoneinfo->start_point_R[i], 
                         gridzoneinfo->start_point_Z[i]);
  }
}

// void update_GridZone_pol_norm_distrb(GridZoneInfo* gridzoneinfo, const double* norm_distb, int n_point)
// {
//   if (!gridzoneinfo || !norm_distb || n_point < 1)
//   {
//     fprintf(stderr, "Invalid input to update_GridZone_pol_norm_distrb.\n");
//     exit(EXIT_FAILURE);
//   }
//   if (gridzoneinfo->npoint != -1 || gridzoneinfo->norm_pol_dist)
//   {
//     fprintf(stderr, "GridZoneInfo already contains norm_pol_dist data.\n");
//     exit(EXIT_FAILURE);
//   }
//   gridzoneinfo->norm_pol_dist = malloc(n_point * sizeof(double));
//   if (!gridzoneinfo->norm_pol_dist)
//   {
//     fprintf(stderr, "Memory allocation failed for norm_pol_dist.\n");
//     exit(EXIT_FAILURE);
//   }
//   for (int i = 0; i < n_point; ++i)
//   {
//     gridzoneinfo->norm_pol_dist[i] = norm_distb[i];
//   }
//   gridzoneinfo->npoint = n_point;
// }
// void update_GridZone_first_pol_points(GridZoneInfo* gridzoneinfo, DLListNode* head)
// {
//   if (!gridzoneinfo || !head)
//   {
//     fprintf(stderr, "Invalid input to update_GridZone_end_curve.\n");
//     exit(EXIT_FAILURE);
//   }
//   DLListNode* node = head;
//   int n = 0;
//   while(node)
//   {
//     n++;
//     node=node->next;
//   }
//   gridzoneinfo->first_pol_points = create_curve(n);
//   if (!gridzoneinfo->first_pol_points) 
//   {
//     fprintf(stderr, "Memory allocation failed for first_pol_points.\n");
//     exit(EXIT_FAILURE);
//   }
//   node=head;
//   for(int i=0;i<n;i++)
//   {
//     add_last_point_curve(gridzoneinfo->first_pol_points,node->r,node->z);
//     node=node->next;
//   }
// }

void write_polsegms_from_dgtrg(DivGeoTrg* trg, const char* filename)
{
  if (!trg || !filename) {
    fprintf(stderr, "Invalid input to write_PolSegmsInfo.\n");
    return;
  }
  int n_region = trg->n_region;
  PolSegmsInfo* polseginfo = create_PolSegmsInfo(n_region);
  size_t len = strlen(trg->topo);
  polseginfo->topo=(char*)malloc(len + 1);
  if (!polseginfo->topo) 
  {
    fprintf(stderr, "Failed to allocate memory for topo.\n");
    exit(EXIT_FAILURE);
  }
  strcpy(polseginfo->topo, trg->topo);

  //According to topology, determine the fowlloing
  if(strncmp(polseginfo->topo,"SNL",3)==0&&n_region==3)
  {
    for(int i=0; i< n_region; ++i)
    {
      polseginfo->xptidx[i]=0; //first X-point
      polseginfo->segmidx[i]=0; //only one segment of sep lines for SN.
      polseginfo->reverse_segm[i]=0; //Not reverse.
    }
    polseginfo->seplineidx[0]=1;
    polseginfo->seplineidx[1]=0;
    polseginfo->seplineidx[2]=2;
  }
  else
  {
    //todo for other configurations
    fprintf(stderr, "Currently only support SNL.\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0;i<n_region;i++)
  {
    if (!trg->zones[i] || !trg->zones[i]->norm_dist || !trg->zones[i]->name) 
    {
      fprintf(stderr, "Invalid region at index %d\n", i);
      exit(EXIT_FAILURE);
    }
    int n_points = trg->zones[i]->n_points;
    polseginfo->polsegments[i]=create_PolSegStr(n_points, 
                                                trg->zones[i]->name);
    for(int j=0; j<n_points; j++)
    {
      polseginfo->polsegments[i]->norm_dist[j]= trg->zones[i]->norm_dist[j];
    }
  }
  write_PolSegmsInfo(polseginfo,filename);
  free_PolSegmsInfo(polseginfo);
}