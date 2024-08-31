#include "datastructure.h"

DLListNode* create_DLListNode(double r, double z)
{
  DLListNode* ptr = NULL;
  DLListNode* new_node = malloc(sizeof(DLListNode));
  if (new_node == NULL){
    printf("Can not create new DLListNode, No memory avaible\n");
    exit(1);
  }
  new_node->r = r;
  new_node->z = z;
  new_node->prev = NULL;
  new_node->next = NULL;
  ptr = new_node;
  return ptr;
};

// head is the head of the double linked list
DLListNode* insert_DLList_at_head(DLListNode* head, double r, double z)
{
  DLListNode* new_node = malloc(sizeof(DLListNode));
  if (new_node == NULL){
    printf("Can not insert DLListNode at head, No memory avaible\n");
    exit(1);
  }
  new_node->r = r;
  new_node->z = z;
  new_node->prev = NULL;
  new_node->next = head;

  if( head != NULL)
  {
    head->prev = new_node;
  }
  else
  {
    printf("The head of DLListNode is empty!\n");
    exit(1);
  }
  return new_node;
}

// free a double linked list
void free_DLList(DLListNode* head)
{
 DLListNode* tmp;
 while(head != NULL )
 {
    tmp = head->next;
    free(head);
    head = tmp;
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