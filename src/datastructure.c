#include "datastructure.h"

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

