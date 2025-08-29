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

//get the DLList from the file
DLListNode* load_DLList_from_file(const char* filename);


// head is the head node of the double linked list
void add_DLListnode_at_head(DLListNode** head_ref, double r, double z);

// get the end node of a double linked list, becareful to use, because the DLL will updated!!!
DLListNode* get_DLList_tailnode(DLListNode* head);

// insert a point at the end of a double linked list
void add_DLListnode_at_tail(DLListNode** prt_end, double r, double z);

DLListNode* copy_DLList(DLListNode* head);

//directly connect the DLList head2 to head1. NOT COPY a new one then connection
//if the tail of head1 is same with the head2, skip=1 will skip the head2 and start from head2.next;
void connect_DLList(DLListNode* head1,DLListNode** head2, int skip);

// free a double linked list
void free_DLList(DLListNode* head);

// print a double linked list
void print_DLList(DLListNode* head);

// write a double linked list to a file
void write_DLList(DLListNode* head, const char* title);

//insert a node between a and b, we assume b=a.next.
void insert_between(DLListNode* a, DLListNode* b, double r, double z);  
// check whether two DLList are intersected or not.
// intersec return 0 else 1
int has_intersection_DLList(DLListNode* head1, DLListNode* head2);

// check whether two DLList are intersected or not.
// intersec return 0 else 1
// idx1 and idx2 are the indexs in head1 and head2 to indicate the no. of segments.
// if no intersection, idx1 and idx2 are reset to -1
int has_intersection_DLList_indexs(DLListNode* head1, DLListNode* head2, int* idx1, int* idx2);

//assunme only 1 intersection point
// r and z are the intersection point
int insert_intersections_DLList(DLListNode* head1, DLListNode* head2, double* r, double* z);

// return how many node are delete, cut the NEXT point of (r,z), NOT (r,z) itself
// -1 means not in the line, 0 means r z is the last
int cut_DLList_after_point(DLListNode* head, double r, double z);

// Cut all nodes before the specified point (r, z), update head to new position
// Returns the number of nodes deleted
// If point not found, returns 0 and head remains unchanged
int cut_DLList_before_point(DLListNode** head, double r, double z);

// split the DLList from point r,z and create a new node with rz and bound to new_head;
int split_DLList_at_point(DLListNode* head, double r, double z, DLListNode** new_head);

//Reverse the DLList
void reverse_DLList(DLListNode** head);

//
double total_length_DLList(DLListNode* head);

// return 0: sucessfully insert the point
// else return 1;
//only work for the situation that the point is confirmed in the ddl line
int insert_point_for_DLList(DLListNode* head, double r, double z);


/*
Following are functions that allocate dynamic memory for 2D,3D,4D arrays
*/

// allocate dynamic 2D array and retunre the pointer
double** allocate_2d_array(const int d1, const int d2);
// free 2D array
void free_2d_array(double **array);

//todo
void write_2d_array(const int d1, const int d2, double** array, const char* filename);


// allocate dynamic 2D array and retunre the pointer
double*** allocate_3d_array(const int d1, const int d2, const int d3);
// free 3D array
void free_3d_array(double ***array);

// Write 3D array to file
// dim=1: write in d1-major order (i->j->k)
// dim=2: write in d2-major order (j->i->k)
void write_3d_array(const int d1, const int d2, const int d3, double*** array, 
                    const char* filename, int dim);
#endif