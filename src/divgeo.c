#include "divgeo.h"
#include "datastructure.h"
#include <stdlib.h>
#include <math.h>
#include "target.h"
#include <stdbool.h>
#define MAX_ITER_DG 1000
#define EPSILON_DG 1.0E-10

// DGClosedStruc* create_DGClosedStruc(int n_point)
// {
//   DGClosedStruc* struc=malloc(sizeof(DGClosedStruc));
//   if(!struc)
//   {
//     fprintf(stderr, "Error: failed to allocate memory for DGClosedStruc\n");
//     return NULL;
//   }
//   struc->n_point=n_point;
//   struc->points=allocate_2d_array(n_point,2);
//   if (!struc->points)
//   {
//     fprintf(stderr, "Error: failed to allocate memory for points in DGClosedStruc\n");
//     free(struc);
//     return NULL;
//   }
//   return struc;
// }

void free_DGClosedStruc(DGClosedStruc* struc)
{
  if (!struc) return;
  free_2d_array(struc->points);
  free(struc);
}

DGRegion* create_DGRegion(int n_level, const char* name)
{
  DGRegion* region = malloc(sizeof(DGRegion));
  if (!region) {
    fprintf(stderr, "Error: failed to allocate memory for DGRegion\n");
    return NULL;
  }

  // Allocate memory for name
  region->name = malloc(strlen(name) + 1);  // +1 for null terminator
  if (!region->name) 
  {
    fprintf(stderr, "Error: failed to allocate memory for region name\n");
    free(region);
    return NULL;
  }
  strcpy(region->name, name);  // Now safe


  region->n_level = n_level;
  region->level = malloc(sizeof(double) * n_level);
  if (!region->level) {
    fprintf(stderr, "Error: failed to allocate memory for region levels\n");
    free(region);
    return NULL;
  }

  return region;
}

void free_DGRegion(DGRegion* region)
{
  if (!region) return;

  if (region->name) 
  {
    free(region->name);
    region->name = NULL;
  }

  if (region->level) 
  {
    free(region->level);
    region->level = NULL;
  }

  free(region);
}

DGZone* create_DGZone(int n_points, const char* name)
{
  if(!name)
  {
    fprintf(stderr, "Error: The name is NULL\n");
    return NULL;
  }
  DGZone* zone = malloc(sizeof(DGZone));
  if (!zone) 
  {
    fprintf(stderr, "Error: failed to allocate memory for DGZone\n");
    return NULL;
  }

  zone->name = malloc(strlen(name) + 1);
  if (!zone->name) 
  {
    fprintf(stderr, "Error: failed to allocate memory for zone name\n");
    free(zone);
    return NULL;
  }
  strcpy(zone->name, name);

  zone->n_points = n_points;
  zone->norm_dist = malloc(sizeof(double) * n_points);
  if (!zone->n_points) 
  {
    fprintf(stderr, "Error: failed to allocate memory for zone levels\n");
    free(zone->name);
    free(zone);
    return NULL;
  }

  return zone;
}

void free_DGZone(DGZone* zone)
{
  if (!zone) return;

  if (zone->name) {
    free(zone->name);
    zone->name = NULL;
  }

  if (zone->norm_dist) {
    free(zone->norm_dist);
    zone->norm_dist = NULL;
  }
  free(zone);
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

static DGRegion* read_region(FILE* fp, int n_level, const char* name)
{
    char line[256];
    DGRegion* region = create_DGRegion(n_level, name);
    region->level[0]=NAN;
    fgets(line, sizeof(line), fp);//skip the line which is 'level'
    for(int i=1; i<n_level; i++)
    {
      fgets(line, sizeof(line), fp);
      line[strcspn(line, "\r\n")] = '\0';
      sscanf(line, "%lf",&(region->level[i]));
    }
    printf("%s\n",region->name);
    printf("n_level: %d\n",n_level);
    for(int i=0; i<n_level; i++)
    {
      printf("%lf\n", region->level[i]);
    }
    return region;
}
static DGZone* read_zone(FILE* fp, int n_points, const char* name)
{
    char line[256];

    DGZone* zone = create_DGZone(n_points, name);
    zone->norm_dist[0]=0.00;
    zone->norm_dist[n_points-1]=1.00;

    fgets(line, sizeof(line), fp);//skip the line which is 'points'

    for(int i=1; i<n_points-1; i++)
    {
      fgets(line, sizeof(line), fp);
      line[strcspn(line, "\r\n")] = '\0';
      sscanf(line, "%lf", &(zone->norm_dist[i]));
    }
    printf("%s", zone->name);
    printf("n_points: %d\n",n_points);

    for(int i=0; i<n_points; i++)
    {
      printf("%lf\n", zone->norm_dist[i]);
    }
    return zone;
}

static Curve* read_curve(FILE* fp, int* n_curve_point)
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

    Curve* curve = create_curve(n_point);
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
    trg->target_curves = count_targets > 0 ? malloc(count_targets * sizeof(Curve*)) : NULL;
    trg->regions = count_regions > 0 ? malloc(count_regions * sizeof(DGRegion*)) : NULL;
    trg->zones = count_zones > 0 ? malloc(count_zones * sizeof(DGZone*)) : NULL;

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
        free_curve(trg->target_curves[i]);
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
        free_DGRegion(trg->regions[i]);
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
        free_DGZone(trg->zones[i]);
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
//DLListNode* head is TargetDDListCurve head
static void cal_points_from_psi(double *psi,double *r, double *z, int n,
                             DLListNode* head,
                             Equilibrium* equ,
                             interpl_2D_1f interp2d1f)
{
  double psi_start, psi_end;
  interp2d1f(r[0],z[0], equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_start, NULL, NULL, NULL);
  interp2d1f(r[n-1],z[n-1], equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_end, NULL, NULL, NULL);

  double psi_head, psi_tail;
  DLListNode* end = get_DLList_endnode(head);
  interp2d1f(head->r,head->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_head, NULL, NULL, NULL);
  interp2d1f(end->r,end->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_tail, NULL, NULL, NULL);

  if((psi_start>psi_head-EPSILON_DG&&psi_start<psi_tail+EPSILON_DG) 
      ||(psi_start<psi_head+EPSILON_DG&&psi_start>psi_tail-EPSILON_DG))
  {
    printf("DEBUG psi_start %lf is in the range psi_head %lf psi_tail %lf\n", psi_start, psi_head, psi_tail);
  }
  else
  {
    fprintf(stderr, "The psi_start %lf is out of range psi_head %lf psi_tail %lf\n!", psi_start, psi_head, psi_tail);
  }

  if( (psi_end>psi_head-EPSILON_DG&&psi_end<psi_tail+EPSILON_DG) 
      ||(psi_end<psi_head+EPSILON_DG&&psi_end>psi_tail-EPSILON_DG) )
  {
    printf("DEBUG psi_start %lf is in the range psi_head %lf psi_tail %lf\n", psi_end, psi_head, psi_tail);
  }
  else
  {
    fprintf(stderr, "The psi_start %lf is out of range psi_head %lf psi_tail %lf\n!", psi_end, psi_head, psi_tail);
  }


  for(int i=0; i<n; i++)
  {
    double psi_prev;
    double psi_curr;
    double psi_next;
    DLListNode* current = head;
    while(current->next!=NULL)
    {
      interp2d1f(head->r,head->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_prev, NULL, NULL, NULL);
      interp2d1f(head->next->r,head->next->z, equ->nw, equ->r, equ->nh, equ->z, equ->psi, &psi_curr, NULL, NULL, NULL);
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
      if (fabs(denom) < 1e-14) 
      {
        printf("Warning: function difference too small in Secant iteration.\n");
        break;
      }
      t_next = t_curr - (psi_curr - psi_target) * (t_curr - t_prev) / denom;
      double r_next = head->r + (head->next->r-head->r)*t_next;
      double z_next = head->z + (head->next->z-head->z)*t_next;
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



void write_dgtrg_to_sn_input(DivGeoTrg* trg, Equilibrium* equ, SeparatrixStr* sep, GradPsiLineStr* gradpsilines)
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
*    STEP1 create TargetDDListCurve
***********************************************/
  int idx_inner_tgt = 1; //start from 0, the second target is the inner target

  TargetDLListCurve* inner_tgt_curve=create_target_curve_from_dgtrg(trg, idx_inner_tgt);
  printf("Inner Target name: %s\n", inner_tgt_curve->name);

/***********************************************
*    STEP2 Sort index for sep and gradpsi lines
***********************************************/
  sort_sep_gradpsiline_by_targetcurve(inner_tgt_curve, sep, gradpsilines);


/*************************************************
*    STEP3 Create the target curve for SOL and PFR
*************************************************/
  double itsct_r, itsct_z;
  // printf("DEBUG sep->index[0] %d\n",sep->index[0]);
  // write_DDList(sep->line_list[sep->index[0]],"debug_sep1");
  // write_DDList(inner_tgt_curve->head,"debug_it");

  insert_intersections_DDList(sep->line_list[sep->index[0]],
                              inner_tgt_curve->head,
                              &itsct_r, &itsct_z);
  
  cut_intersections_DDList(sep->line_list[sep->index[0]],itsct_r, itsct_z);
  // printf("DEBUG intersect:%d\n", 
  //       has_intersection_DDList(sep->line_list[sep->index[0]],
  //                               inner_tgt_curve->head));
  
  TargetDLListCurve* sol_tgt_curve=create_target_curve();
  split_intersections_target_curve(inner_tgt_curve, itsct_r, itsct_z, sol_tgt_curve);
  reverse_DDList_in_target_curve(inner_tgt_curve);
  //need to update the number and name!
  change_name_target_curve(inner_tgt_curve, "PFR");
  change_name_target_curve(sol_tgt_curve, "SOL");

  printf("%s target curve has %d points\n", inner_tgt_curve->name, inner_tgt_curve->n);
  printf("%s target curve has %d points\n", sol_tgt_curve->name, sol_tgt_curve->n);

  write_DDList(sep->line_list[sep->index[0]],"debug_sep1");
  write_DDList(inner_tgt_curve->head,"debug_pfr");
  write_DDList(sol_tgt_curve->head,"debug_sol");

/*************************************************
*    STEP4 Calculate the start points used for trcing
*************************************************/



/*************************************************
*    Free Memory 
*************************************************/
  free_target_curve(inner_tgt_curve);
  free_target_curve(sol_tgt_curve);

}


