#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

// head is the head node of the double linked list
void insert_DLList_at_head(DLListNode** prt_head, double r, double z);

// get the end node of a double linked list, becareful to use, because the DLL will updated!!!
DLListNode* end_DLListNode(DLListNode* head);

// insert a point at the end of a double linked list
void insert_DLList_at_end(DLListNode** prt_end, double r, double z);

// free a double linked list
void free_DLList(DLListNode* head);

// print a double linked list
void print_DLList(DLListNode* head);

// write a double linked list to a file
void write_DDList(DLListNode* head, const char* title);

#endif