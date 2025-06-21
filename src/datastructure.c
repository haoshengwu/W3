#include "datastructure.h"
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#define EPS_DATASTR 1.0E-10

//************** Double linked list related ******************//
DLListNode* create_DLListNode(double r, double z)
{
  
  DLListNode* new_node = malloc(sizeof(DLListNode));
  if (new_node == NULL){
    fprintf(stderr, "Cannot allocate DLListNode: Out of memory.\n");
    exit(EXIT_FAILURE);
    return NULL;
  }
  new_node->r = r;
  new_node->z = z;
  new_node->prev = NULL;
  new_node->next = NULL;
  return new_node;
};

// head is the head of the double linked list
void add_DLListnode_at_head(DLListNode** head_ref, double r, double z)
{
  DLListNode* new_node = create_DLListNode(r, z);
  if (new_node == NULL)
  {
    exit(EXIT_FAILURE);
  }

  if (*head_ref == NULL)
  {
      new_node->next = NULL;
  }
  else
  {
      new_node->next = (*head_ref);
      (*head_ref)->prev = new_node;
  }
  *head_ref = new_node;
}

// get the end node of a double linked list, becareful to use, because the DLL will updated!!!
DLListNode* get_DLList_tailnode(DLListNode* head)
{
  if (head == NULL)
  {
    fprintf(stderr, "Empty input for get_DLList_tailnode.\n");
    exit(EXIT_FAILURE);
  }

  while (head->next != NULL)
  {
    head = head->next;
  }
  return head;
}

// insert a point r,z to the end of double linked list and return new end node
void add_DLListnode_at_tail(DLListNode** end_ref, double r, double z)
{
  if (end_ref == NULL || *end_ref == NULL) 
  {
    fprintf(stderr, "end_ref or *end_ref is NULL!\n");
    exit(EXIT_FAILURE);
  }
  
  DLListNode* new_node = create_DLListNode(r, z);
  if (!new_node) exit(EXIT_FAILURE);

  new_node->prev = *end_ref;
  (*end_ref)->next = new_node;
  *end_ref = new_node;

}

DLListNode* copy_DLList(DLListNode* head) 
{
  if (head == NULL) return NULL;

  DLListNode* new_head = NULL;
  DLListNode* new_tail = NULL;

  DLListNode* current = head;
  while (current != NULL) {
    DLListNode* new_node = (DLListNode*)malloc(sizeof(DLListNode));
    if (!new_node) 
    {
      fprintf(stderr,"malloc failed");
      exit(EXIT_FAILURE);
    }

    new_node->r = current->r;
    new_node->z = current->z;
    new_node->next = NULL;
    new_node->prev = new_tail;

    if (new_tail != NULL) 
    {
      new_tail->next = new_node;
    } 
    else 
    {
      new_head = new_node;  
    }

    new_tail = new_node;
    current = current->next;
  }
    return new_head;
}

void connect_DLList(DLListNode* head1,DLListNode** head2, int skip)
{
  if(!head1||!head2||!(*head2))
  {
    fprintf(stderr,"Empty input for connect_DLList.\n");
    exit(EXIT_FAILURE);
  }
  DLListNode* tail1=get_DLList_tailnode(head1);
  DLListNode* cur = *head2;
  
  int delete_duplicate=0;
  if(fabs(tail1->r-cur->r)<EPS_DATASTR&&
      fabs(tail1->z-cur->z)<EPS_DATASTR)
      {
        if(skip==1)
        {
          cur=cur->next;
          delete_duplicate=1;
          printf("Duplicate point %lf %lf and will be delted.\n", tail1->r, tail1->z);
        }
        else
        {
          printf("WARNING: the head2 R Z is same with tail1 R Z.\n");
          printf("WARNING: they are two same points in the DLList\n");
        }
      }
  if (cur) 
  {
    cur->prev = tail1;
    tail1->next = cur;
  }

  if(delete_duplicate)
  {
    (*head2)->next=NULL;
    (*head2)->prev=NULL;
    free(*head2);
    *head2=NULL;
  }
}


// free a double linked list
void free_DLList(DLListNode* head)
{
  DLListNode* next;
  for (DLListNode* tmp = head; tmp != NULL; tmp = next)
  {
    next = tmp->next;
    free(tmp);
  }
}
// print a double linked list
void print_DLList(DLListNode* head)
{
 while(head != NULL )
 {
    printf("%lf, %lf\n", head->r, head->z);
    head = head->next;
 }
 printf("Already print all value in the double linked list");
}

// write a double linked list to a file
void write_DLList(DLListNode* head, const char* filename)
{
    FILE* file = fopen(filename, "w");
    if (!file) 
    {
      fprintf(stderr, "Can not open %s\n", filename);
      return;
    }
    
    while (head !=NULL) 
    {
        // the unit is mm, it can be changed in the future
        fprintf(file, "%.12f  %.12f\n", head->r, head->z);  
        head = head->next;
    }

    fclose(file);
    printf("write the values in %s\n", filename);
}

/* 2D cross product: (B - A) × (C - A) */
static inline double cross(double ax, double ay, 
                           double bx, double by,
                           double cx, double cy) 
{
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

/* Check if (cx, cy) lies on segment (ax, ay)-(bx, by) */
static inline int on_segment(double ax, double ay, 
                             double bx, double by, 
                             double cx, double cy) 
{
    return fmin(ax, bx) - EPS_DATASTR <= cx && cx <= fmax(ax, bx) + EPS_DATASTR &&
           fmin(ay, by) - EPS_DATASTR <= cy && cy <= fmax(ay, by) + EPS_DATASTR;
}

/* Return 0 if segments intersect (including endpoints), 1 if not */
static int has_intersection_segment(double x1, double y1, double x2, double y2,
                             double x3, double y3, double x4, double y4)
{
    // Step 1: Bounding box check
    if (fmax(x1, x2) + EPS_DATASTR < fmin(x3, x4) ||
        fmax(x3, x4) + EPS_DATASTR < fmin(x1, x2) ||
        fmax(y1, y2) + EPS_DATASTR < fmin(y3, y4) ||
        fmax(y3, y4) + EPS_DATASTR < fmin(y1, y2)) 
    {
        return 1;
    }

    // Step 2: Cross product signs
    double d1 = cross(x3, y3, x4, y4, x1, y1);
    double d2 = cross(x3, y3, x4, y4, x2, y2);
    double d3 = cross(x1, y1, x2, y2, x3, y3);
    double d4 = cross(x1, y1, x2, y2, x4, y4);

    if ((d1 * d2 < 0) && (d3 * d4 < 0)) return 0;

    if (fabs(d1) < EPS_DATASTR && on_segment(x3, y3, x4, y4, x1, y1)) return 0;
    if (fabs(d2) < EPS_DATASTR && on_segment(x3, y3, x4, y4, x2, y2)) return 0;
    if (fabs(d3) < EPS_DATASTR && on_segment(x1, y1, x2, y2, x3, y3)) return 0;
    if (fabs(d4) < EPS_DATASTR && on_segment(x1, y1, x2, y2, x4, y4)) return 0;

    return 1;
}

// 0 means have intersection, 1 means no.
int has_intersection_DLList(DLListNode* head1, DLListNode* head2) 
{
  for (DLListNode* a = head1; a && a->next; a = a->next) 
    {
      for (DLListNode* b = head2; b && b->next; b = b->next)
         {
            if (has_intersection_segment(
                    a->r, a->z, a->next->r, a->next->z,
                    b->r, b->z, b->next->r, b->next->z)==0) 
            {
                return 0;  // Found
            }
        }
    }
    return 1; // no intersection
}

// Compute intersection point of two segments (if any)
// Only works if the intersection is a single point (not overlapping segments)
// intersect return 0, otherwise return 1


// Determine whether two line segments intersect and calculate the intersection point
static int get_intersection(double x1, double y1, double x2, double y2,
                           double x3, double y3, double x4, double y4,
                           double* intsect_x, double* intsect_y)
{
    // First check if they intersect
    // Bounding box rejection
    if (fmax(x1, x2) + EPS_DATASTR < fmin(x3, x4) ||
        fmax(x3, x4) + EPS_DATASTR < fmin(x1, x2) ||
        fmax(y1, y2) + EPS_DATASTR < fmin(y3, y4) ||
        fmax(y3, y4) + EPS_DATASTR < fmin(y1, y2)) 
    {
        *intsect_x=NAN;
        *intsect_y=NAN;
        return 1;
    }

    // Cross products
    double d1 = cross(x3, y3, x4, y4, x1, y1);
    double d2 = cross(x3, y3, x4, y4, x2, y2);
    double d3 = cross(x1, y1, x2, y2, x3, y3);
    double d4 = cross(x1, y1, x2, y2, x4, y4);

    // General case
    if ((d1 * d2 < 0) && (d3 * d4 < 0)) {
        // Use parametric intersection: compute t
        double dx1 = x2 - x1;
        double dy1 = y2 - y1;
        double dx2 = x4 - x3;
        double dy2 = y4 - y3;

        double denom = dx1 * dy2 - dy1 * dx2;
        if (fabs(denom) < EPS_DATASTR) 
        {
            *intsect_x=NAN;
            *intsect_y=NAN;
            return 1; // Parallel or degenerate
        }
        double t = ((x3 - x1) * dy2 - (y3 - y1) * dx2) / denom;

        *intsect_x = x1 + t * dx1;
        *intsect_y = y1 + t * dy1;
        return 0;
    }

    // Special cases: endpoints overlap
    if (fabs(d1) < EPS_DATASTR && fmin(x3,x4)-EPS_DATASTR <= x1 && x1 <= fmax(x3,x4)+EPS_DATASTR &&
                                 fmin(y3,y4)-EPS_DATASTR <= y1 && y1 <= fmax(y3,y4)+EPS_DATASTR) 
    {
        *intsect_x = x1; *intsect_y = y1; return 0;
    }
    if (fabs(d2) < EPS_DATASTR && fmin(x3,x4)-EPS_DATASTR <= x2 && x2 <= fmax(x3,x4)+EPS_DATASTR &&
                                 fmin(y3,y4)-EPS_DATASTR <= y2 && y2 <= fmax(y3,y4)+EPS_DATASTR) 
    {
        *intsect_x = x2; *intsect_y = y2; return 0;
    }
    if (fabs(d3) < EPS_DATASTR && fmin(x1,x2)-EPS_DATASTR <= x3 && x3 <= fmax(x1,x2)+EPS_DATASTR &&
                                 fmin(y1,y2)-EPS_DATASTR <= y3 && y3 <= fmax(y1,y2)+EPS_DATASTR) 
    {
        *intsect_x = x3; *intsect_y = y3; return 0;
    }
    if (fabs(d4) < EPS_DATASTR && fmin(x1,x2)-EPS_DATASTR <= x4 && x4 <= fmax(x1,x2)+EPS_DATASTR &&
                                 fmin(y1,y2)-EPS_DATASTR <= y4 && y4 <= fmax(y1,y2)+EPS_DATASTR) 
    {
        *intsect_x = x4; *intsect_y = y4; return 0;
    }
    return 1;
}






void insert_between(DLListNode* a, DLListNode* b, double r, double z) 
{
    if (!a || !b) 
    {
        fprintf(stderr, "Empty input for insert_between.\n");
        exit(EXIT_FAILURE);
    }

    if (a->next != b || b->prev != a)
    {
        fprintf(stderr, "Non-contiguous nodes in insert_between.\n");
        exit(EXIT_FAILURE);
    }

    DLListNode* new_node = create_DLListNode(r, z);
    if (!new_node) return;

    new_node->prev = a;
    new_node->next = b;
    a->next = new_node;
    b->prev = new_node;
}

// Returns 0 if an intersection is found and inserted, 1 otherwise
int insert_intersections_DLList(DLListNode* head1, DLListNode* head2, double* r_ptr, double* z_ptr) 
{
    int found = 1;
    for (DLListNode* a = head1; a && a->next; a = a->next) 
    {
      for (DLListNode* b = head2; b && b->next; b = b->next) 
      {
        if (get_intersection(
                  a->r, a->z, a->next->r, a->next->z,
                  b->r, b->z, b->next->r, b->next->z,
                  r_ptr, z_ptr)==0) 
        {
          // Insert the intersection point into both DLLists
          // Use two separate nodes to avoid sharing the same memory
          if( (fabs(a->r-*r_ptr)<EPS_DATASTR && fabs(a->z-*z_ptr)<EPS_DATASTR)||
              (fabs(a->next->r-*r_ptr)<EPS_DATASTR && fabs(a->next->z-*z_ptr)<EPS_DATASTR))
          {
            printf("DEBUG: %.12f %.12f already in the DLList1, not need to insert\n",*r_ptr, *z_ptr);
          }
          else
          {
            insert_between(a, a->next, *r_ptr, *z_ptr);  // insert into list 1
          }
          if( (fabs(b->r-*r_ptr)<EPS_DATASTR && fabs(b->z-*z_ptr)<EPS_DATASTR)||
              (fabs(b->next->r-*r_ptr)<EPS_DATASTR && fabs(b->next->z-*z_ptr)<EPS_DATASTR))
          {
            printf("DEBUG: %.12f %.12f already in the DLList2, not need to insert\n",*r_ptr, *z_ptr);
          }
          else
          {
            insert_between(b, b->next, *r_ptr, *z_ptr);  // insert into list 2
          }
          found=0;
          break; // prevent multiple insertions on the same segment
        }
      }
      if(found==0) 
      {
        printf("DEBUG intersection point: %.12f %.12f\n", *r_ptr, *z_ptr);
        break; // Break outer loop after successful insertion
      }
    }
    if(found)
    {
      printf("WARNING: No intersection point found.\n");
    }
    return found;
}

int cut_DLList_from_intersections(DLListNode* head, double r, double z) {
    DLListNode* current = head;
    bool found = false;
    // Step 1: Find the node with matching (r, z)
    while (current) {
    if (fabs(current->r - r) < EPS_DATASTR && fabs(current->z - z) < EPS_DATASTR) 
        {
          break;
        }
        found = true;
        current = current->next;
    }

    if (!found) {
        printf("The point (%.12f, %.12f) is not in the DLList.\n", r, z);
        return 0;
    }

    // Step 2: Delete all nodes after current
    DLListNode* node = current->next;
    current->next = NULL; // Cut the list here
    if (node) node->prev = NULL; // Disconnect the next node (optional safety)

    int count = 0;
    while (node) {
        DLListNode* next = node->next;
        free(node);
        node = next;
        count++;
    }
    printf("Cut %d points in the DLList.\n", count);
    return count; // Number of nodes deleted
}


// Split the list at the first (r,z) point — clone (r,z) node as new head
// 0 means slit
int split_intersections_DLList(DLListNode* head, double r, double z, DLListNode** new_head) 
{
    DLListNode* curr = head;
    while (curr) {
        if (fabs(curr->r - r) < EPS_DATASTR && fabs(curr->z - z) < EPS_DATASTR) {
            DLListNode* cloned = create_DLListNode(r, z);
            if (!cloned) return 0;

            // Set up new_head
            *new_head = cloned;

            // Connect cloned node to the rest of original list
            cloned->next = curr->next;
            if (cloned->next) {
                cloned->next->prev = cloned;
            }

            // Cut off the original list at curr
            curr->next = NULL;

            return 0;
        }
        curr = curr->next;
    }

    *new_head = NULL;
    return 1;
}

void reverse_DLList(DLListNode** head)
{
    // Check if the list is NULL or empty
    if (!head || !(*head)) return;

    DLListNode* current = *head;  // Current node in traversal
    DLListNode* prev = NULL;      // Previous node (will become next after reversal)
    DLListNode* next = NULL;      // Temporarily holds the original next node

    // Traverse the list and reverse the next and prev pointers
    while (current != NULL)
    {
        next = current->next;     // Save original next node

        current->next = prev;     // Reverse the next pointer
        current->prev = next;     // Reverse the prev pointer

        prev = current;           // Move prev to current
        current = next;           // Move to next node in original order
    }

    // After the loop, prev is the new head
    *head = prev;
}


double total_length_DLList(DLListNode* head)
{
  if (!head || !head->next)
  {
    fprintf(stderr, "Empty or too short input for total_length_DLList.\n");
    exit(EXIT_FAILURE);
  }

  double length = 0.0;
  DLListNode* cur = head;

  while (cur->next)
  {
    length += hypot(cur->r - cur->next->r, cur->z - cur->next->z);
    cur = cur->next;
  }

  return length;
}

//*****************************************************************


double **allocate_2d_array(const int d1, const int d2)
{
/*
allocate dynamic 2D array and retunre the pointer
d1: first dimension number
d2: second dimenion number
*/
  if (d1 <= 0 || d2 <= 0)
  {
    fprintf(stderr, "Error: Invalid dimensions for 2D array!\n");
    return NULL;
  }

  double **array = (double **)malloc(d1 * sizeof(double *));
  if (array == NULL)
  {
    fprintf(stderr, "Error: Memory allocation failed for pointer array!\n");
    return NULL;
  }

  array[0] = (double *)malloc(d2*d1*sizeof(double));
  if (array[0] == NULL)
  {
    fprintf(stderr, "Error: Memory allocation failed for pointer array[0]!\n");
    return NULL;
  }

  for (int i = 0; i<d1; i++)
  {
    array[i] = array[0] + i*d2;
  }
  
  return array;
}

void free_2d_array(double **array)
{
  free(array[0]);
  free(array);
} 



double*** allocate_3d_array(const int d1, const int d2, const int d3)
{
/*
allocate dynamic 2D array and retunre the pointer
d1: first dimension number
d2: second dimension number
d3: third dimension number
*/
  if (d1 <= 0 || d2 <= 0 || d3 <= 0)
  {
    fprintf(stderr, "Error: Invalid dimensions!\n");
    return NULL;
  }

  double ***array = (double ***)malloc(d1 * sizeof(double **));
  if (array == NULL)
  {
    fprintf(stderr, "Error: Memory allocation failed for pointer array!\n");
    return NULL;
  }

  array[0] = (double **)malloc(d2 * d1 * sizeof(double *));
  if (array[0] == NULL)
  {
    fprintf(stderr, "Error: Memory allocation failed for pointer array array[0]!\n");
    free(array);
    return NULL;
  }

  array[0][0] = (double *)malloc(d3 * d2 * d1 * sizeof(double));
  if (array[0][0] == NULL)
  {
    fprintf(stderr, "Error: Memory allocation failed for pointer array array[0][0]!\n");
    free(array[0]); // Free the second dimension pointers
    free(array);
    return NULL;
  }

  for(int i=0; i<d1; i++)
  {
    array[i] = array[0] + i*d2;
    for(int j=0; j<d2; j++)
    {
      array[i][j] = array[0][0] + (i*d2 + j)*d3;
    }
  }

  return array;
}

void free_3d_array(double ***array)
{
  free(array[0][0]);
  free(array[0]);
  free(array);
}

