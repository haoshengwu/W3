#include "target.h"

TargetCurve* create_target_curve()
{
  TargetCurve* tgt_cur = malloc(sizeof(TargetCurve));
  memset(tgt_cur->name, 0, sizeof(tgt_cur->name));  
  tgt_cur->n=-1;
  tgt_cur->head = NULL;
  return tgt_cur;
}

int add_point_target_curve(TargetCurve* tgt_cur, double r, double z)
{
  if(!tgt_cur)
  {
    fprintf(stderr,"Empty of TargetCurve!\n");
    exit(EXIT_FAILURE);
  }

  if(!tgt_cur->head)
  {
    tgt_cur->head=create_DLListNode(r,z);
    tgt_cur->n=0;
    return 0;
  }
  else
  {
    DLListNode* tail = get_DLList_endnode(tgt_cur->head);
    insert_DLList_at_end(&tail, r, z);
    tgt_cur->n++;
    return 0;
  }
}
int change_name_target_curve(TargetCurve* curve, const char* name)
{
  if(!curve || !name)
  {
    fprintf(stderr,"Wrong input for change_name_target_curve!\n");
    return 1;
  }
  strncpy(curve->name, name, sizeof(curve->name) - 1);
  curve->name[sizeof(curve->name) - 1] = '\0';
  return 0;
}

TargetCurve* create_target_curve_from_dgtrg(DivGeoTrg* trg, int n)
{
  TargetCurve* tgt_cur=create_target_curve();
  if (!tgt_cur) 
  {
    fprintf(stderr, "Failed to allocate memmory for target_curve");
    exit(EXIT_FAILURE);
  }
  char name[32]="Target";
  snprintf(name, sizeof(name), "TargetCurve_%d\n", n);
  change_name_target_curve(tgt_cur, name);
  int n_point=trg->n_target_curve[n];
  //  printf("Target Curve Name: %s\n", tgt_cur->name);
  for(int i=0; i<n_point;i++)
  {
    // printf("DEBUG r=%.4f z=%.4f\n", r, z);
    double r = trg->target_curves[n]->points[i][0];
    double z = trg->target_curves[n]->points[i][1];
    // printf("DEBUG i=%d\n", i);
    printf("DEBUG r=%.4f z=%.4f\n", r, z);
    add_point_target_curve(tgt_cur, r, z);
  }
  return tgt_cur;
}

void free_target_curve(TargetCurve* curve)
{
  if(!curve->head)
  {
    return;
  }
  free_DLList(curve->head);
  free(curve);
}

TargetCurve* reverse_target_curve(TargetCurve** curve)
{
  if (!curve || !(*curve) || !(*curve)->head) {
    fprintf(stderr, "Input curve is NULL or empty in reverse_and_free_target_curve!\n");
    return NULL;
  }

  TargetCurve* original = *curve;

  TargetCurve* reversed = create_target_curve();
  change_name_target_curve(reversed, original->name);


  DLListNode* tail = get_DLList_endnode(original->head);
  while (tail) {
    add_point_target_curve(reversed, tail->r, tail->z);
    tail = tail->prev;
  }

  free_target_curve(original);
  *curve = reversed;

  return reversed;
}


void printf_target_curve(TargetCurve* curve)
{
  if (!curve) return;

  printf("Target Curve Name: %s", curve->name);
  printf("The points R Z:\n");

  DLListNode* node = curve->head;
  while (node)
  {
    printf("%lf %lf\n", node->r, node->z);
    node = node->next;
  }
}