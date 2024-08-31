#ifndef __DATASTRUCTURE_H
#define __DATASTRUCTURE_H

#include <stdlib.h>
#include <stdio.h>
// similar to DG ngroups, there are some simple data structure are defined here.

// a double-linked list, in each node, there is the coordinate of a point (r,z)

typedef struct DLList
{
  double r,z;
  struct DLList* prev;
  struct DLList* next;
}
DLListNode;

DLListNode* create_DLListNode(double r, double z);

// head is the head of the double linked list
DLListNode* insert_DLList_at_head(DLListNode* head, double r, double z);

// free a double linked list
void free_DLList(DLListNode* head);

// print a double linked list
void print_DLList(DLListNode* head);

#endif