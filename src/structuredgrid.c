#include "structuredgrid.h"
#include "datastructure.h"
#include <stdlib.h>



GridZone* allocate_GridZone()
{
  GridZone* z=malloc(sizeof(GridZone));
  if (!z) {
    fprintf(stderr, "Error: failed to allocate memory for GridZone\n");
    return NULL;
  }
  z->name=NULL;
  z->np=-1;
  z->nr=-1;
  
  z->grid_curveset=NULL;

  z->start_point_R=NULL;
  z->start_point_Z=NULL;
  z->guard_start=NULL;
  z->guard_end=NULL;
  z->pasmin = NULL;
  z->norm_pol_dist=NULL;
 
  z->first_boundary=NULL;
  z->second_boundary=NULL;

  z->target_curve=NULL;

  return z;
}

int load_GridZone_from_file(GridZone* z, const char* filename)
{
  if (!z || !filename) 
  {
    fprintf(stderr, "Invalid arguments to load_GridZone_from_file\n");
    return 1;
  }
  
  // Save the filename into z->name
  z->name = malloc(strlen(filename) + 1);
  if (!z->name) {
    fprintf(stderr, "Memory allocation failed for GridZone name\n");
    return 1;
  }
  strcpy(z->name, filename);

  FILE* fp = fopen(filename, "r");
  if (!fp) 
  {
    fprintf(stderr, "Fail to open the file: %s\n", filename);
    goto fail;
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
    goto fail;
  }

  // Allocate memory
  z->start_point_R = malloc(sizeof(double) * z->nr);
  z->start_point_Z = malloc(sizeof(double) * z->nr);

  z->guard_start = malloc(sizeof(double) * z->nr);
  z->guard_end = malloc(sizeof(double) * z->nr);
  z->pasmin = malloc(sizeof(double) * z->nr);
  z->norm_pol_dist = malloc(sizeof(double) * z->np);

  if (!z->start_point_R || !z->start_point_Z
      || !z->guard_start || !z->guard_end || !z->pasmin || !z->norm_pol_dist) 
  {
    fprintf(stderr, "Memory allocation failed for start_points or guard arrays.\n");
    goto fail;
  }

  // Skip to next valid line for distribution
  while (fgets(line, sizeof(line), fp))
  {
    if (line[0] == '#' || line[0] == '\n') continue;
    break;
  }

  for (int i = 0; i<z->np; i++)
  {
    if (sscanf(line, "%lf", &z->norm_pol_dist[i]) != 1) 
    {
      fprintf(stderr, "Failed to read norm_pol_dist[%d] from file.\n", i);
      goto fail;
    }
    if (!fgets(line, sizeof(line), fp)) 
    {
      fprintf(stderr, "Unexpected EOF while reading norm_pol_dist.\n");
      goto fail;
    }
  }

// Skip to next valid line for distribution
  while (line[0] == '#' || line[0] == '\n') 
  {
    if (!fgets(line, sizeof(line), fp)) 
    {
      fprintf(stderr, "Unexpected EOF before tracing data.\n");
      goto fail;
    }
}

  // Read data for each tracing line
  for (int i = 0; i < z->nr; i++) 
  {
    if (sscanf(line, "%lf %lf %lf %lf %lf",
              &(z->start_point_R[i]), &(z->start_point_Z[i]),
              &(z->guard_start[i]), &(z->guard_end[i]), &(z->pasmin[i])) != 5)
    {
      fprintf(stderr,
              "Invalid format at line %d for tracing input (expecting 5 values).\n", i + 1);
      goto fail;
    }
    if (!fgets(line, sizeof(line), fp) && i < z->nr - 1) 
    {
      fprintf(stderr, "Unexpected EOF when reading tracing points.\n");
      goto fail;
    }
  }

  if (!feof(fp)) 
  {
    printf("Warning: File %s has additional lines after expected input.\n", filename);
  }
  fclose(fp);
  return 0;
  
fail:
    if (fp) fclose(fp);
    if (z->name) { free(z->name); z->name = NULL; }
    if (z->start_point_R) { free(z->start_point_R); z->start_point_R = NULL; }
    if (z->start_point_Z) { free(z->start_point_Z); z->start_point_Z = NULL; }
    if (z->guard_start) { free(z->guard_start); z->guard_start = NULL;}
    if (z->guard_end) {free(z->guard_end);  z->guard_end = NULL;}
    if (z->pasmin) {free(z->pasmin);     z->pasmin = NULL;}
    if (z->norm_pol_dist) {free(z->norm_pol_dist); z->norm_pol_dist = NULL;}
    return 1;
}

int load_GridZone_first_boundary(GridZone* z, Curve* first_boundary)
{
    if (!z || !first_boundary) 
    {
        fprintf(stderr, "Invalid input to load_GridZone_first_boundary.\n");
        return 1;
    }

    if (z->first_boundary == NULL)
        z->first_boundary = create_curve(first_boundary->n_point);

    return copy_curve(z->first_boundary, first_boundary);
}

int load_GridZone_second_boundary(GridZone* z, Curve* second_boundary)
{
    if (!z || !second_boundary) 
    {
        fprintf(stderr, "Invalid input to load_GridZone_second_boundary.\n");
        return 1;
    }

    if (z->second_boundary == NULL)
        z->second_boundary = create_curve(second_boundary->n_point);

    return copy_curve(z->second_boundary, second_boundary);
}

int load_GridZone_target_curve(GridZone* z, Curve* target_curve)
{
    if (!z || !target_curve) 
    {
        fprintf(stderr, "Invalid input to load_GridZone_target_curve.\n");
        return 1;
    }

    if (z->target_curve == NULL)
        z->target_curve = create_curve(target_curve->n_point);

    return copy_curve(z->target_curve, target_curve);
}

void free_GridZone(GridZone** z)
{
  if (!z) return;

  //free name
  free((*z)->name);
  (*z)->name=NULL;

  free_curveset((*z)->grid_curveset);
  (*z)->grid_curveset=NULL;

  free((*z)->start_point_R);
  (*z)->start_point_R=NULL;

  free((*z)->start_point_Z);
  (*z)->start_point_Z=NULL;

  free((*z)->guard_start);
  (*z)->guard_start=NULL;

  free((*z)->guard_end);
  (*z)->guard_end=NULL;

  free((*z)->pasmin);
  (*z)->pasmin=NULL;

  free((*z)->norm_pol_dist);
  (*z)->norm_pol_dist=NULL;

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

