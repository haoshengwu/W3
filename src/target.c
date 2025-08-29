#include "target.h"

TargetDLListCurve* create_target_curve()
{
  TargetDLListCurve* tgt_cur = malloc(sizeof(TargetDLListCurve));
  memset(tgt_cur->name, 0, sizeof(tgt_cur->name));  
  tgt_cur->n=-1;
  tgt_cur->head = NULL;
  return tgt_cur;
}

int add_point_target_curve(TargetDLListCurve* tgt_cur, double r, double z)
{
  if(!tgt_cur)
  {
    fprintf(stderr,"Empty of TargetDLListCurve!\n");
    exit(EXIT_FAILURE);
  }

  if(!tgt_cur->head)
  {
    tgt_cur->head=create_DLListNode(r,z);
    tgt_cur->n=1;
    return 0;
  }
  else
  {
    DLListNode* tail = get_DLList_tailnode(tgt_cur->head);
    add_DLListnode_at_tail(&tail, r, z);
    tgt_cur->n++;
    return 0;
  }
}

int change_name_target_curve(TargetDLListCurve* tgt_cur, const char* name)
{
  if(!tgt_cur || !name)
  {
    fprintf(stderr,"Wrong input for change_name_target_curve!\n");
    return 1;
  }
  strncpy(tgt_cur->name, name, sizeof(tgt_cur->name) - 1);
  tgt_cur->name[sizeof(tgt_cur->name) - 1] = '\0';
  return 0;
}


void free_target_curve(TargetDLListCurve* tgt_cur)
{
  if(!tgt_cur->head)
  {
    return;
  }
  free_DLList(tgt_cur->head);
  free(tgt_cur);
}


void reverse_DLList_in_target_curve(TargetDLListCurve* tgt_cur)
{
  if (!tgt_cur || !(tgt_cur->head)) 
  {
    fprintf(stderr, "Input curve is NULL or empty in reverse_and_free_target_curve!\n");
    return;
  }
  reverse_DLList(&tgt_cur->head);
}

void split_intersections_target_curve(TargetDLListCurve* tgt_cur, 
                                      double r, double z, 
                                      TargetDLListCurve* new_tgt_cur)
{
  int status;
  status=split_DLList_at_point(tgt_cur->head, r, z, &(new_tgt_cur->head));
  if(status)
  {
    fprintf(stderr,"Failed to split the DLList in the target curve.\n");
    exit(EXIT_FAILURE);
  }
  update_number_target_curve(tgt_cur);
  update_number_target_curve(new_tgt_cur);
}


void cut_target_curve(TargetDLListCurve* tgt_cur, double r, double z)
{
  if(!tgt_cur||!tgt_cur->head)
  {
    fprintf(stderr, "Empty input of cut_target_curve");
    exit(EXIT_FAILURE);
  }
  
  int status;
  status = cut_DLList_after_point(tgt_cur->head, r ,z);
  
  if(status==0)
  {
    fprintf(stderr,"The point is not in the target curve!\n");
    exit(EXIT_FAILURE);
  }
}



void update_number_target_curve(TargetDLListCurve* tgt_cur)
{
  if(!tgt_cur)
  {
    fprintf(stderr,"Empty input for update_number_target_curve.\n");
    exit(EXIT_FAILURE);
  }
  int n=0;
  DLListNode* head=tgt_cur->head;
  while(head!=NULL)
    {
      n=n+1;
      head=head->next;
    }
  tgt_cur->n=n;
}

void printf_target_curve(TargetDLListCurve* tgt_cur)
{
  if (!tgt_cur) return;

  printf("Target Curve Name: %s", tgt_cur->name);
  printf("The points R Z:\n");

  DLListNode* node = tgt_cur->head;
  while (node)
  {
    printf("%lf %lf\n", node->r, node->z);
    node = node->next;
  }
}

void sort_sep_gradpsiline_by_targetcurve(TargetDLListCurve* inner_tgt_cur,
                                         TargetDLListCurve* outer_tgt_cur,
                                         SeparatrixStr* sep,
                                         GradPsiLineStr* gradpsilines)
{
  if(!inner_tgt_cur->head || !outer_tgt_cur->head)
  {
    fprintf(stderr,"Target Curve is NULL!\n");
    exit(EXIT_FAILURE);
  }

  DLListNode* inner_head = inner_tgt_cur->head;
  DLListNode* outer_head = outer_tgt_cur->head;

  int start=-1;
  for(int i=0; i<4; i++)
  {
    int sep_idx1=-1;
    int sep_idx2=-1;
    int tmp1=-1;
    int tmp2=-1;
    
    has_intersection_DLList_indexs(sep->line_list[i], inner_head, &sep_idx1, &tmp1);
    has_intersection_DLList_indexs(sep->line_list[i], outer_head, &sep_idx2, &tmp2);

    // printf("DEBUG i sep_idx1 sep_idx2 %d %d %d\n", i, sep_idx1, sep_idx2);

    if(sep_idx1>=0 && tmp1>=0) // There is intersection between the sep line and inner target
    {
      if(sep_idx2==-1 && tmp2==-1) //The sep line do not intersect with outer target
      {
        start=i;
      }
      else if(sep_idx2>=0 && tmp2 >=0) //The sep line also intersect with outer target
      {
        if(sep_idx1<sep_idx2)
        {
          start=i;
          break;
        }
      }
    }
  }

  if(start==-1)
  {
    fprintf(stderr, "No intersection between target curve and separatrix line!\n");
    exit(EXIT_FAILURE);
  }

  for(int i=0; i<4; i++)
  {
    sep->index[i]=(start+i)%4;
    gradpsilines->index[i] = (start+i)%4;
  }

#ifdef DEBUG
  for(int i=0; i<4; i++)
  {
    printf("Separatrix Line %d, Index: %d\n", i, sep->index[i]);
  }

  for(int i=0; i<4; i++)
  {
    printf("Grad PSI Line %d, Index: %d\n", i, gradpsilines->index[i]);
  }
#endif

}

TargetDLListCurve* create_core_curve_from_gradpsilines(GradPsiLineStr* gradpsilines, int idx)
{
  if (!gradpsilines || !gradpsilines->line_list || !gradpsilines->index) {
    fprintf(stderr, "Invalid gradpsilines or its members are NULL\n");
    exit(EXIT_FAILURE);
  }

  if (idx < 0 || idx >= 4) {
    fprintf(stderr, "Index %d is out of bounds (nlines = %d)\n", idx, 4);
    exit(EXIT_FAILURE);
  }

  TargetDLListCurve* core_curve = create_target_curve();
  if (!core_curve) {
    fprintf(stderr, "Failed to allocate memory for target_curve\n");
    exit(EXIT_FAILURE);
  }

  char name[32];
  snprintf(name, sizeof(name), "Core");
  change_name_target_curve(core_curve, name);

  DLListNode* head = gradpsilines->line_list[gradpsilines->index[idx]];
  core_curve->head = copy_DLList(head);

  if (!core_curve->head) {
    fprintf(stderr, "Failed to copy DLList for core curve\n");
    exit(EXIT_FAILURE);
  }

  update_number_target_curve(core_curve);
  return core_curve;
}