#include "datastructure.h"
#include <math.h>
#include <stdio.h>
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
void insert_DLList_at_head(DLListNode** head_ref, double r, double z)
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
DLListNode* get_DLList_endnode(DLListNode* head)
{
  if (head == NULL)
  {
    fprintf(stderr, "The double linked list is empty!\n");
    return NULL;
  }

  while (head->next != NULL)
  {
    head = head->next;
  }
  return head;
}

// insert a point r,z to the end of double linked list and return new end node
void insert_DLList_at_end(DLListNode** end_ref, double r, double z)
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
void write_DDList(DLListNode* head, const char* filename)
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
        fprintf(file, "%lf  %lf\n", head->r, head->z);  
        head = head->next;
    }

    fclose(file);
    printf("write the values in %s\n", filename);
}

// Determines the orientation of the ordered triplet (x1, y1), (x2, y2), (x3, y3)
// Returns:
//   0 → Colinear (all points lie on a straight line)
//   1 → Clockwise turn (right turn)
//   2 → Counterclockwise turn (left turn)
static int orientation(double x1, double y1, double x2, double y2, double x3, double y3) 
{
    double val = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);

    if (val > 0) return 1;  // Clockwise
    if (val < 0) return 2;  // Counterclockwise
    return 0;               // Colinear
}


// Checks if point (x2, y2) lies on the line segment from (x1, y1) to (x3, y3)
// Assumes the three points are colinear
static int on_segment(double x1, double y1, double x2, double y2, double x3, double y3) 
{
  if (x2 <= fmax(x1, x3) && x2 >= fmin(x1, x3) &&
      y2 <= fmax(y1, y3) && y2 >= fmin(y1, y3))
      {
        return 0;
      }
      return 1;
}

// Returns 0 if the two segments (x1, y1)-(x2, y2) and (x3, y3)-(x4, y4) intersect
// Returns 1 otherwise
static int segments_intersect(double x1, double y1, double x2, double y2,
                       double x3, double y3, double x4, double y4) {
    // Find orientations of the four combinations
    int o1 = orientation(x1, y1, x2, y2, x3, y3);
    int o2 = orientation(x1, y1, x2, y2, x4, y4);
    int o3 = orientation(x3, y3, x4, y4, x1, y1);
    int o4 = orientation(x3, y3, x4, y4, x2, y2);

    // General case: two endpoints of one segment lie on opposite sides of the other segment
    if (o1 != o2 && o3 != o4)
        return 0;

    // Special colinear cases: check if endpoints lie on the other segment
    if (o1 == 0 && on_segment(x1, y1, x3, y3, x2, y2)==0) return 0;
    if (o2 == 0 && on_segment(x1, y1, x4, y4, x2, y2)==0) return 0;
    if (o3 == 0 && on_segment(x3, y3, x1, y1, x4, y4)==0) return 0;
    if (o4 == 0 && on_segment(x3, y3, x2, y2, x4, y4)==0) return 0;

    // No intersection found
    return 1;
}

// 0 means have intersection, 1 means no.
int has_intersection_DDList(DLListNode* head1, DLListNode* head2) 
{
  for (DLListNode* a = head1; a && a->next; a = a->next) 
    {
      for (DLListNode* b = head2; b && b->next; b = b->next)
         {
            if (segments_intersect(
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
static int compute_intersection_point(
    double x1, double y1, double x2, double y2,
    double x3, double y3, double x4, double y4,
    double* ix, double* iy)
{
    // Calculate denominators
    double denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);

    if (fabs(denom) < 1e-10) {
        // Lines are parallel or coincident
        return 1;
    }

    // Compute intersection using determinant formula
    double px_num = (x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4);
    double py_num = (x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4);

    double px = px_num / denom;
    double py = py_num / denom;

    // Check if the intersection point lies on both segments
    if (px < fmin(x1,x2) - 1e-10 || px > fmax(x1,x2) + 1e-10 ||
        px < fmin(x3,x4) - 1e-10 || px > fmax(x3,x4) + 1e-10 ||
        py < fmin(y1,y2) - 1e-10 || py > fmax(y1,y2) + 1e-10 ||
        py < fmin(y3,y4) - 1e-10 || py > fmax(y3,y4) + 1e-10) {
        return 1; // Intersection outside the segments
    }

    *ix = px;
    *iy = py;
    return 0;
}

void insert_between(DLListNode* a, DLListNode* b, double r, double z) 
{
    DLListNode* new_node = create_DLListNode(r, z);
    if (!new_node) return;
    new_node->prev = a;
    new_node->next = b;
    if (a) a->next = new_node;
    if (b) b->prev = new_node;
}


//0 means intersect, else 1 means not.
int insert_intersections_DDList(DLListNode* head1, DLListNode* head2, double* r_ptr, double* z_ptr) 
{
    int intersect = 1;
    for (DLListNode* a = head1; a && a->next; a = a->next) {
        for (DLListNode* b = head2; b && b->next; b = b->next) {
              if (compute_intersection_point(
                    a->r, a->z, a->next->r, a->next->z,
                    b->r, b->z, b->next->r, b->next->z,
                    r_ptr, z_ptr)==0) {
                // Insert the intersection point into both DLLists
                // Use two separate nodes to avoid sharing the same memory
                insert_between(a, a->next, *r_ptr, *z_ptr);  // insert into list 1
                insert_between(b, b->next, *r_ptr, *z_ptr);  // insert into list 2
                intersect=0;
                break; // prevent multiple insertions on the same segment
            }
        }
      if(intersect==0) break;
    }
    return intersect; // total number of intersections inserted
}


  
int cut_intersections_DDList(DLListNode* head, double r, double z) {
    DLListNode* current = head;
    const double tol = 1e-10;

    // Step 1: Find the node with matching (r, z)
    while (current) {
        if (fabs(current->r - r) < tol && fabs(current->z - z) < tol) {
            break;
        }
        current = current->next;
    }

    if (!current) {
        fprintf(stderr, "The point (%.10f, %.10f) is not in the DLList.\n", r, z);
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

    return count; // Number of nodes deleted
}

// Split the list at the first (r,z) point — clone (r,z) node as new head
int split_intersections_DDList(DLListNode* head, double r, double z, DLListNode** new_head) 
{
    const double tol = 1e-10;
    DLListNode* curr = head;

    while (curr) {
        if (fabs(curr->r - r) < tol && fabs(curr->z - z) < tol) {
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

            return 1;
        }
        curr = curr->next;
    }

    *new_head = NULL;
    return 0;
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

