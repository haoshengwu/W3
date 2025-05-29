#include "structuredgrid.h"
#include "datastructure.h"
#include <stdlib.h>

static int copy_curve(Curve* c1, Curve* c2)
{
    if (!c1 || !c2) {
        fprintf(stderr, "Null pointer passed to copy_curve.\n");
        return 1;
    }

    if (c1->n_point != c2->n_point) {
        fprintf(stderr, "The size of two curves is not the same!\n");
        return 1;
    }

    for (size_t i = 0; i < c1->n_point; i++) {
        c1->points[i][0] = c2->points[i][0];
        c1->points[i][1] = c2->points[i][1];
    }
    return 0;
}

Curve* create_curve(size_t n_point)
{   
  Curve* c = malloc(sizeof(Curve));
  if (!c)
  {
    fprintf(stderr, "Error: Memory allocation failed for a curve!\n");
    return NULL;
  } 
  c->n_point = n_point;
  c->points = allocate_2d_array(n_point,2);
  return c;
}

void free_curve(Curve* c) {
    if (!c) return;
    free_2d_array(c->points);  // Free row pointers
    free(c);         // Free the Curve struct
}

CurveSet* create_curveset(size_t n_curve, size_t n_point) {
    CurveSet* cs = malloc(sizeof(CurveSet));
    if (!cs) {
        fprintf(stderr, "Error: Failed to allocate CurveSet\n");
        return NULL;
    }

    cs->n_curve = n_curve;
    cs->curves = malloc(sizeof(Curve*) * n_curve);

    for (size_t i = 0; i < n_curve; i++) {
        cs->curves[i] = create_curve(n_point);  // Return Curve*
        if (!cs->curves[i]) {
            fprintf(stderr, "Error: Failed to create curve %zu\n", i);
            // Free previously allocated curves before exiting
            for (size_t j = 0; j < i; j++)
                free_curve(cs->curves[j]);
            free(cs->curves);
            free(cs);
            return NULL;
        }
    }

    return cs;
}

void free_curveset(CurveSet *cs)
{
  if (!cs) return;
  for(size_t i=0; i<cs->n_curve;i++)
  {
    free_curve(cs->curves[i]);
  }
  free(cs->curves);
  free(cs);
}


Zone* allocate_zone()
{
  Zone* z=malloc(sizeof(Zone));
  if (!z) {
    fprintf(stderr, "Error: failed to allocate memory for Zone\n");
    return NULL;
  }
  z->name=NULL;
  z->np=-1;
  z->nr=-1;
  
  z->zone_grid=NULL;

  z->start_points=NULL;
  z->guard_head=NULL;
  z->guard_end=NULL;
  z->pasmin = NULL;
  z->distribution[0]=0.00;
  z->distribution[1]=0.00;

  z->first_boundary=NULL;
  z->second_boundary=NULL;

  z->target_curve=NULL;

  return z;
}

int load_zone_from_file(Zone* z, const char* filename)
{
  if (!z || !filename) 
  {
    fprintf(stderr, "Invalid arguments to load_zone_from_file\n");
    return 1;
  }
  
  // Save the filename into z->name
  z->name = malloc(strlen(filename) + 1);
  if (!z->name) {
    fprintf(stderr, "Memory allocation failed for zone name\n");
    return 1;
  }
  strcpy(z->name, filename);

  FILE* fp = fopen(filename, "r");
  if (!fp) 
  {
    fprintf(stderr, "Fail to open the file: %s\n", filename);
    free(z->name);
    z->name = NULL;
    return 1;
  }

  char line[256];

  while (fgets(line, sizeof(line), fp))
  {
    if (line[0] == '#' || line[0] == '\n') continue;
    break;
  }

  // Read np and nr
  if (sscanf(line, "%d %d", &(z->np), &(z->nr)) != 2) 
  {
    fprintf(stderr, "Invalid format for np and nr in file: %s\n", filename);
    fclose(fp);
    free(z->name);
    z->name = NULL;
    return 1;
  }

  // Allocate memory
  z->start_points = allocate_2d_array(z->nr, 2);
  z->guard_head = malloc(sizeof(double) * z->nr);
  z->guard_end = malloc(sizeof(double) * z->nr);
  z->pasmin = malloc(sizeof(double) * z->nr);

  if (!z->start_points || !z->guard_head || !z->guard_end || !z->pasmin) 
  {
    fprintf(stderr, "Memory allocation failed for start_points or guard arrays.\n");
    fclose(fp);
    return 1;
  }

  // Skip to next valid line for distribution
  while (fgets(line, sizeof(line), fp))
  {
    if (line[0] == '#' || line[0] == '\n') continue;
    break;
  }

  if (sscanf(line, "%lf %lf", &(z->distribution[0]), &(z->distribution[1])) != 2) 
  {
    fprintf(stderr, "Invalid format for distribution values.\n");
    fclose(fp);
    return 1;
  }
  printf("deltp1: %lf, deltpn: %lf\n", z->distribution[0], z->distribution[1]);

  // Read data for each tracing line
  for (int i = 0; i < z->nr; i++) 
  {
    if (!fgets(line, sizeof(line), fp)) 
    {
      fprintf(stderr, "Unexpected end of file while reading start_point at index %d.\n", i);
      fclose(fp);
      return 1;
    }
    if (sscanf(line, "%lf %lf %lf %lf %lf",
              &(z->start_points[i][0]), &(z->start_points[i][1]),
              &(z->guard_head[i]), &(z->guard_end[i]), &(z->pasmin[i])) != 5)
        {
          fprintf(stderr,
                "Invalid format at line %d for tracing input (expecting 5 values).\n", i + 1);
          fclose(fp);
          return 1;
        }
    }
  if (!feof(fp)) 
  {
    printf("Warning: File %s has additional lines after expected input.\n", filename);
  }
  fclose(fp);
  return 0;
}

int load_zone_first_boundary(Zone* z, Curve* first_boundary)
{
    if (!z || !first_boundary) 
    {
        fprintf(stderr, "Invalid input to load_zone_first_boundary.\n");
        return 1;
    }

    if (z->first_boundary == NULL)
        z->first_boundary = create_curve(first_boundary->n_point);

    return copy_curve(z->first_boundary, first_boundary);
}

int load_zone_second_boundary(Zone* z, Curve* second_boundary)
{
    if (!z || !second_boundary) 
    {
        fprintf(stderr, "Invalid input to load_zone_second_boundary.\n");
        return 1;
    }

    if (z->second_boundary == NULL)
        z->second_boundary = create_curve(second_boundary->n_point);

    return copy_curve(z->second_boundary, second_boundary);
}

int load_zone_target_curve(Zone* z, Curve* target_curve)
{
    if (!z || !target_curve) 
    {
        fprintf(stderr, "Invalid input to load_zone_target_curve.\n");
        return 1;
    }

    if (z->target_curve == NULL)
        z->target_curve = create_curve(target_curve->n_point);

    return copy_curve(z->target_curve, target_curve);
}

void free_zone(Zone** z)
{
  if (!z) return;

  //free name
  free((*z)->name);
  (*z)->name=NULL;

  free_curveset((*z)->zone_grid);
  (*z)->zone_grid=NULL;

  free_2d_array((*z)->start_points);
  (*z)->start_points=NULL;

  free((*z)->guard_head);
  (*z)->guard_head=NULL;

  free((*z)->guard_end);
  (*z)->guard_end=NULL;

  free((*z)->pasmin);
  (*z)->pasmin=NULL;

  free_curve((*z)->first_boundary);
  (*z)->first_boundary=NULL;

//NOT always has senond boudanry curve;
  if((*z)->second_boundary)
  {
    free_curve((*z)->second_boundary);
    (*z)->second_boundary=NULL;
  }

  free_curve((*z)->target_curve);
  (*z)->target_curve=NULL;

  free(*z); 
  *z = NULL;
}