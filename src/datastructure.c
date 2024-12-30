#include "datastructure.h"

//************** Double linked list related ******************//
DLListNode* create_DLListNode(double r, double z)
{
  
  DLListNode* new_node = malloc(sizeof(DLListNode));
  if (new_node == NULL){
    printf("Can not create new DLListNode, No memory avaible\n");
    return NULL;
  }
  new_node->r = r;
  new_node->z = z;
  new_node->prev = NULL;
  new_node->next = NULL;
  return new_node;
};

// head is the head of the double linked list
void insert_DLList_at_head(DLListNode** ptr_head, double r, double z)
{
  DLListNode* new_node = malloc(sizeof(DLListNode));
  if (new_node == NULL){
    printf("Can not insert DLListNode at head, No memory avaible\n");
    exit(1);
  }
    new_node->r = r;
    new_node->z = z;
    new_node->prev = NULL;

  if (*ptr_head == NULL)
  {
      new_node->next = NULL;
  }
  else
  {
      new_node->next = (*ptr_head);
      (*ptr_head)->prev = new_node;
  }
  *ptr_head = new_node;
}

// get the end node of a double linked list, becareful to use, because the DLL will updated!!!
DLListNode* end_DLListNode(DLListNode* head)
{
  if (head == NULL)
  {
    printf("the Double linked list is emptry!\n");
    return NULL;
  }

  while (head->next != NULL)
  {
    head = head->next;
  }
  return head;
}

// insert a point r,z to the end of double linked list and return new end node
void insert_DLList_at_end(DLListNode** ptr_end, double r, double z)
{
  if (ptr_end == NULL || (*ptr_end) == NULL)
  {
    printf("the ptr_end or end node is emptry!\n");
    exit(1);
  }

  DLListNode* new_node = malloc(sizeof(DLListNode));
  
  if (new_node == NULL){
    printf("Can not insert DLListNode at end, No memory avaible\n");
    exit(1);
  }
    new_node->r = r;
    new_node->z = z;
    new_node->next = NULL;
    new_node->prev = (*ptr_end);
    (*ptr_end)->next = new_node;
    (*ptr_end) = new_node;
}

// free a double linked list
void free_DLList(DLListNode* head)
{
 DLListNode* tmp = head;
 DLListNode* next;
 while(tmp != NULL )
 {
    next =tmp->next;
    free(tmp);
    tmp = next;
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
    fprintf(stderr, "Error: Invalid dimensions!\n");
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

