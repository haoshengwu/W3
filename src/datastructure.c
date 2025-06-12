#include "datastructure.h"
#include <math.h>
#include <stdio.h>
#define EPSILON_DATASTR 1.0E-10

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

// Compute the cross product of two vectors
static double cross(double x1, double y1, double x2, double y2) {
  return x1 * y2 - x2 * y1;
}

// Check if a value is within [0, 1] considering floating-point tolerance EPSILON
static int inRange(double value) {
    return value >= -EPSILON_DATASTR && value <= 1.0 + EPSILON_DATASTR;
}

// Returns 0 if the two segments (x1, y1)-(x2, y2) and (x3, y3)-(x4, y4) intersect
// Returns 1 otherwise
static int segments_intersect(double x1, double y1, double x2, double y2,
                       double x3, double y3, double x4, double y4) {
    double dx1 = x2 - x1, dy1 = y2 - y1; // Direction vector of the first segment
    double dx2 = x4 - x3, dy2 = y4 - y3; // Direction vector of the second segment
    double delta = cross(dx1, dy1, dx2, dy2); // Cross product to check for parallelism

    if (fabs(delta) < EPSILON_DATASTR) return 1; // Segments are parallel or colinear

    // Calculate the parameters of intersection point on both segments
    double s = cross(x3 - x1, y3 - y1, dx2, dy2) / delta;
    double t = cross(x3 - x1, y3 - y1, dx1, dy1) / delta;

    // Check if the intersection point lies within both segments
    if (inRange(s) && inRange(t)) 
    {
        return 0; // Segments intersect
     }
    // printf("WARNING: NO INTERSECTION FOUND.\n");
    return 1; // Segments do not intersect
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


// Determine whether two line segments intersect and calculate the intersection point
static int getIntersection(double x1, double y1, double x2, double y2,
                    double x3, double y3, double x4, double y4,
                    double *ix, double *iy) {
    double dx1 = x2 - x1, dy1 = y2 - y1; // Direction vector of the first segment
    double dx2 = x4 - x3, dy2 = y4 - y3; // Direction vector of the second segment
    double delta = cross(dx1, dy1, dx2, dy2); // Cross product to check for parallelism

    if (fabs(delta) < EPSILON_DATASTR) return 1; // Segments are parallel or colinear

    // Calculate the parameters of intersection point on both segments
    double s = cross(x3 - x1, y3 - y1, dx2, dy2) / delta;
    double t = cross(x3 - x1, y3 - y1, dx1, dy1) / delta;

    // Check if the intersection point lies within both segments
    if (inRange(s) && inRange(t)) 
    {
       *ix = x1 + s * dx1;
       *iy = y1 + s * dy1;
        return 0; // Segments intersect
     }
    // printf("WARNING: NO INTERSECTION FOUND.\n");
    return 1; // Segments do not intersect
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

// Returns 0 if an intersection is found and inserted, 1 otherwise
int insert_intersections_DDList(DLListNode* head1, DLListNode* head2, double* r_ptr, double* z_ptr) 
{
    int found = 1;
    for (DLListNode* a = head1; a && a->next; a = a->next) 
    {
      for (DLListNode* b = head2; b && b->next; b = b->next) 
      {
        if (getIntersection(
                  a->r, a->z, a->next->r, a->next->z,
                  b->r, b->z, b->next->r, b->next->z,
                  r_ptr, z_ptr)==0) 
        {
              // Insert the intersection point into both DLLists
              // Use two separate nodes to avoid sharing the same memory
              insert_between(a, a->next, *r_ptr, *z_ptr);  // insert into list 1
              insert_between(b, b->next, *r_ptr, *z_ptr);  // insert into list 2
              found=0;
              break; // prevent multiple insertions on the same segment
        }
      }
      if(found==0) 
      {
        printf("DEBUG intersection point: %lf %lf\n", *r_ptr, *z_ptr);
        break;
      }
    }
    if(found)
    {
      printf("WARNING: No intersection point found.\n");
    }
    return found;
}

int cut_intersections_DDList(DLListNode* head, double r, double z) {
    DLListNode* current = head;

    // Step 1: Find the node with matching (r, z)
    while (current) {
        if (fabs(current->r - r) < EPSILON_DATASTR && fabs(current->z - z) < EPSILON_DATASTR) {
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
// 0 means slit
int split_intersections_DDList(DLListNode* head, double r, double z, DLListNode** new_head) 
{
    DLListNode* curr = head;
    while (curr) {
        if (fabs(curr->r - r) < EPSILON_DATASTR && fabs(curr->z - z) < EPSILON_DATASTR) {
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

