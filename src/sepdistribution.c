#include "sepdistribution.h"

EdgeSegment* create_EdgeSegment() 
{
  EdgeSegment* segment = malloc(sizeof(EdgeSegment));
  if (!segment) {
    fprintf(stderr, "Failed to allocate memory for EdgeSegment.\n");
    exit(EXIT_FAILURE);
  }

  segment->tail = NULL;
  segment->head = NULL;
  segment->n_size = -1;
  segment->norm_dist = NULL;
  segment->gridpoint_curve = NULL;
  segment->next = NULL;
  return segment;
}

SepDistStr* create_SepDistStr_from_sep(SeparatrixStr* sep)
{
  SepDistStr* sepdist=malloc(sizeof(SepDistStr));
    if (!sepdist) {
    fprintf(stderr, "Failed to allocate memory for SepDistStr.\n");
    exit(EXIT_FAILURE);
  }
  sepdist->xpt_r=sep->xpt_r;
  sepdist->xpt_z=sep->xpt_z;
  for(int i=0; i<4; i++)
  {
    sepdist->index[i]=-1;
    sepdist->edges[i]=create_EdgeSegment();
    sepdist->edges[i]->head = copy_DLList(sep->line_list[i]);
  }
  return sepdist;
}

void free_EdgeSegment(EdgeSegment *seg) 
{
  while (seg) {
    EdgeSegment *next = seg->next;

    if (seg->norm_dist) 
    {
      free(seg->norm_dist);
      seg->norm_dist = NULL;
    }

    if (seg->gridpoint_curve) 
    {
      free_curve(seg->gridpoint_curve);
      seg->gridpoint_curve=NULL;
    }

    if (seg->head) {
      free_DLList(seg->head);
      seg->head = NULL;
      seg->tail = NULL;  
    }
    free(seg);
    seg = next;
  }
}

void free_SepDistStr(SepDistStr* sepdist)
{
  for(int i=0; i<4; i++)
  {
    if(sepdist->edges[i])
    {
      free_EdgeSegment(sepdist->edges[i]);
    }
  }
  free(sepdist);
}


void update_sn_SepDistStr_from_GridZoneInfo(SepDistStr* sepdist, GridZoneInfo* gzinfo)
{
  //===== 1. sort separatrix lines in sepdist
  DLListNode* head_inner=create_DLListNode(gzinfo->start_point_R[0],
                                           gzinfo->start_point_Z[0]);
  DLListNode* tail_inner=head_inner;
  for(int i=1; i<gzinfo->nr; i++)
  {
    add_DLListnode_at_tail(&tail_inner, gzinfo->start_point_R[i],gzinfo->start_point_Z[i]);
  }
  write_DLList(head_inner, "DEBUG_gzinfo_inner");

  int start=-1;
  for(int i=0; i<4; i++)
  {
    if(has_intersection_DLList(sepdist->edges[i]->head, head_inner)==0)
    {
      printf("DEBUG intersection line: %d\n", i);
      start = i;
      break;
    }
  }
  if(start==-1)
  {
    fprintf(stderr, "No intersection between DLLlist and sepdist line!\n");
    exit(EXIT_FAILURE);
  }
  for(int i=0; i<4; i++)
  {
    sepdist->index[i]=(start+i)%4;
    printf("sep index %d: %d\n", i, sepdist->index[i]);
  }

  //===== 2. Cutting sep lines;
  DLListNode* head_outer=create_DLListNode(gzinfo->start_point_R[0],
                                           gzinfo->start_point_Z[0]);
  DLListNode* tail_outer=head_outer;
  for(int i=1; i<gzinfo->end_curve->n_point; i++)
  
  {
    add_DLListnode_at_tail(&tail_outer, 
                          gzinfo->end_curve->points[i].x,
                          gzinfo->end_curve->points[i].y);
  }
  write_DLList(head_outer, "DEBUG_gzinfo_outer");

  double itsct_r_inner, itsct_z_inner;
  int idx = sepdist->index[0];
  insert_intersections_DLList(sepdist->edges[idx]->head,
                              head_inner,
                              &itsct_r_inner, &itsct_z_inner);
  cut_DLList_from_intersections(sepdist->edges[idx]->head,itsct_r_inner, itsct_z_inner);

  double itsct_r_outer, itsct_z_outer;
  idx = sepdist->index[1];
  insert_intersections_DLList(sepdist->edges[idx]->head,
                              head_outer,
                              &itsct_r_outer, &itsct_z_outer);
  cut_DLList_from_intersections(sepdist->edges[idx]->head,itsct_r_outer, itsct_z_outer);

  //===== 3. Free
  free_DLList(head_inner);
  free_DLList(head_outer);
  printf("Finish update SepDistStr by GridZoneInfo\n");
}

void update_sn_SepDistStr_from_PolSegmsInfo(SepDistStr* sepdist, PolSegmsInfo* polseginfo)
{
  if (!sepdist || !polseginfo) {
    fprintf(stderr, "Empty input for update_sn_SepDistStr_from_PolSegmsInfo.\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < 3; ++i) {
    if (!polseginfo->polsegments[i]) {
      fprintf(stderr, "Empty polsegment[%d] for update_sn_SepDistStr_from_PolSegmsInfo.\n", i);
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < 3; ++i) {
    int idx = sepdist->index[i];

    if (!sepdist->edges[idx]) {
      fprintf(stderr, "Empty EdgeSegment[%d] for update_sn_SepDistStr_from_PolSegmsInfo.\n", idx);
      exit(EXIT_FAILURE);
    }

    if (sepdist->edges[idx]->norm_dist != NULL) {
      fprintf(stderr, "norm_dist for edge[%d] is already allocated. Aborting to avoid overwrite.\n", idx);
      exit(EXIT_FAILURE);
    }

    int n_size = polseginfo->polsegments[i]->n_points;
    sepdist->edges[idx]->n_size = n_size;
    sepdist->edges[idx]->norm_dist = malloc(n_size * sizeof(double));
    if (!sepdist->edges[idx]->norm_dist) {
      fprintf(stderr, "Memory allocation failed for norm_dist of edge[%d].\n", idx);
      exit(EXIT_FAILURE);
    }

    for (int j = 0; j < n_size; ++j) {
      sepdist->edges[idx]->norm_dist[j] = polseginfo->polsegments[i]->norm_dist[j];
    }
  }
}