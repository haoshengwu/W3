#include "gridzoneinfo.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>


GridZoneInfo* allocate_GridZoneInfo()
{
  GridZoneInfo* z=malloc(sizeof(GridZoneInfo));
  if (!z) {
    fprintf(stderr, "Error: failed to allocate memory for GridZoneInfo\n");
    return NULL;
  }
  z->topo=NULL;
  z->name=NULL;

  z->nr=-1;
  z->start_point_R=NULL;
  z->start_point_Z=NULL;
  z->guard_start=NULL;
  z->guard_end=NULL;
  z->pasmin = NULL;
  
  z->start_curve=NULL;
  z->end_curve=NULL;

  z->n_polsegm1=-1;
  z->xptidx1=NULL;
  z->seplineidx1=NULL;
  z->segmidx1=NULL;
  z->reverse_segm1=NULL;


  z->n_polsegm2=-1;
  z->xptidx2=NULL;
  z->seplineidx2=NULL;
  z->segmidx2=NULL;
  z->reverse_segm2=NULL;


  return z;
}

void free_GridZoneInfo(GridZoneInfo** z)
{
  if (!z) return;

  //free name

  free((*z)->topo);
  (*z)->topo=NULL;

  free((*z)->name);
  (*z)->name=NULL;

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

  free_curve((*z)->start_curve);
  (*z)->start_curve=NULL;

  free_curve((*z)->end_curve);
  (*z)->end_curve=NULL;

  free((*z)->xptidx1);
  (*z)->xptidx1=NULL;

  free((*z)->seplineidx1);
  (*z)->seplineidx1=NULL;
  
  free((*z)->segmidx1);
  (*z)->segmidx1=NULL;

  free((*z)->reverse_segm1);
  (*z)->reverse_segm1=NULL;

  free((*z)->xptidx2);
  (*z)->xptidx2=NULL;

  free((*z)->seplineidx2);
  (*z)->seplineidx2=NULL;
  
  free((*z)->segmidx1);
  (*z)->segmidx2=NULL;

  free((*z)->reverse_segm2);
  (*z)->reverse_segm2=NULL;

  free(*z); 
  *z = NULL;
}
void write_GridZoneInfo(GridZoneInfo* z, const char* filename)
{
  if (!z || !filename) {
    fprintf(stderr, "Invalid input to write_GridZoneInfo.\n");
    exit(EXIT_FAILURE);
  }

  FILE* fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s for writing.\n", filename);
    exit(EXIT_FAILURE);
  }

  fprintf(fp, "#topo\n%s\n", z->topo ? z->topo : "Untopo");
  fprintf(fp, "#name\n%s\n", z->name ? z->name : "Unnamed");

  fprintf(fp, "#nr\n%d\n", z->nr);
  fprintf(fp, "#start_point_R Start_point_Z guardstart guardend pasmin\n");
  for (int i = 0; i < z->nr; ++i) {
    fprintf(fp, "%.12f %.12f %.12f %.12f %.12f\n",
            z->start_point_R[i], z->start_point_Z[i],
            z->guard_start[i], z->guard_end[i], z->pasmin[i]);
  }

  fprintf(fp, "#start_curve\n");
  if (z->start_curve && z->start_curve->points) {
    fprintf(fp, "%zu\n", z->start_curve->n_point);
    for (size_t i = 0; i < z->start_curve->n_point; ++i) {
      fprintf(fp, "%.12f %.12f\n", z->start_curve->points[i].x, z->start_curve->points[i].y);
    }
  } else {
    fprintf(fp, "0\n# No start_curve defined.\n");
  }

  fprintf(fp, "#end_curve\n");
  if (z->end_curve && z->end_curve->points) {
    fprintf(fp, "%zu\n", z->end_curve->n_point);
    for (size_t i = 0; i < z->end_curve->n_point; ++i) {
      fprintf(fp, "%.12f %.12f\n", z->end_curve->points[i].x, z->end_curve->points[i].y);
    }
  } else {
    fprintf(fp, "0\n# No end_curve defined.\n");
  }

  fprintf(fp, "#xptidx1 seplineidx1 segmidx1 reverse\n");
  fprintf(fp, "%d\n", z->n_polsegm1);
  for (int i = 0; i < z->n_polsegm1; ++i) {
    fprintf(fp, "%d %d %d %d\n", z->xptidx1[i], z->seplineidx1[i], z->segmidx1[i], z->reverse_segm1[i]);
  }

  fprintf(fp, "#xptidx2 seplineidx2 segmidx2 reverse\n");
  if (z->n_polsegm2 < 0) {
    fprintf(fp, "NaN\n");
  } else {
    fprintf(fp, "%d\n", z->n_polsegm2);
    for (int i = 0; i < z->n_polsegm2; ++i) {
      fprintf(fp, "%d %d %d %d\n", z->xptidx2[i], z->seplineidx2[i], z->segmidx2[i], z->reverse_segm2[i]);
    }
  }

  fclose(fp);
  printf("Finish writing GridZoneInfo %s to file %s\n", z->name ? z->name : "Unnamed", filename);
}

GridZoneInfo* load_GridZoneInfo_from_input(const char* filename)
{
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    perror("Failed to open input file");
    exit(EXIT_FAILURE);
  }

  GridZoneInfo* z = allocate_GridZoneInfo();
  char line[256];

  while (fgets(line, sizeof(line), fp)) {
    if (strncmp(line, "#topo", 5) == 0) {
      fgets(line, sizeof(line), fp);
      line[strcspn(line, "\r\n")] = 0;
      size_t len = strlen(line);
      z->topo = malloc(len + 1);
      if (z->topo) strcpy(z->topo, line);
    } else if (strncmp(line, "#name", 5) == 0) {
      fgets(line, sizeof(line), fp);
      line[strcspn(line, "\r\n")] = 0;
      size_t len = strlen(line);
      z->name = malloc(len + 1);
      if (z->name) strcpy(z->name, line);
    } else if (strncmp(line, "#nr", 3) == 0) {
      fgets(line, sizeof(line), fp);
      sscanf(line, "%d", &z->nr);
      int nr = z->nr;
      z->start_point_R = malloc(nr * sizeof(double));
      z->start_point_Z = malloc(nr * sizeof(double));
      z->guard_start = malloc(nr * sizeof(double));
      z->guard_end = malloc(nr * sizeof(double));
      z->pasmin = malloc(nr * sizeof(double));
    } else if (strncmp(line, "#start_point_R", 14) == 0) {
      for (int i = 0; i < z->nr; ++i) {
        fgets(line, sizeof(line), fp);
        sscanf(line, "%lf %lf %lf %lf %lf",
               &z->start_point_R[i], &z->start_point_Z[i],
               &z->guard_start[i], &z->guard_end[i], &z->pasmin[i]);
      }
      } else if (strncmp(line, "#start_curve", 12) == 0) {
      fgets(line, sizeof(line), fp);
      size_t n;
      sscanf(line, "%zu", &n);
      if (n > 0) {
        z->start_curve = create_curve(n);
        for (size_t i = 0; i < n; ++i) {
          fgets(line, sizeof(line), fp);
          double x, y;
          sscanf(line, "%lf %lf", &x, &y);
          add_last_point_curve(z->start_curve, x, y);
        }
      }
    } else if (strncmp(line, "#end_curve", 10) == 0) {
      fgets(line, sizeof(line), fp);
      size_t n;
      sscanf(line, "%zu", &n);
      if (n > 0) {
        z->end_curve = create_curve(n);
        for (size_t i = 0; i < n; ++i) {
          fgets(line, sizeof(line), fp);
          double x, y;
          sscanf(line, "%lf %lf", &x, &y);
          add_last_point_curve(z->end_curve, x, y);
        }
      }
    } else if (strncmp(line, "#xptidx1", 8) == 0) {
      fgets(line, sizeof(line), fp);
      sscanf(line, "%d", &z->n_polsegm1);
      int n = z->n_polsegm1;
      z->xptidx1 = malloc(n * sizeof(int));
      z->seplineidx1 = malloc(n * sizeof(int));
      z->segmidx1 = malloc(n * sizeof(int));
      z->reverse_segm1 = malloc(n * sizeof(int));
      for (int i = 0; i < n; ++i) {
        fgets(line, sizeof(line), fp);
        sscanf(line, "%d %d %d %d",
               &z->xptidx1[i], &z->seplineidx1[i],
               &z->segmidx1[i], &z->reverse_segm1[i]);
      }
    } else if (strncmp(line, "#xptidx2", 8) == 0) {
      fgets(line, sizeof(line), fp);
      if (strncmp(line, "NaN", 3) == 0) {
        z->n_polsegm2 = -1;
      } else {
        sscanf(line, "%d", &z->n_polsegm2);
        int n = z->n_polsegm2;
        z->xptidx2 = malloc(n * sizeof(int));
        z->seplineidx2 = malloc(n * sizeof(int));
        z->segmidx2 = malloc(n * sizeof(int));
        z->reverse_segm2 = malloc(n * sizeof(int));
        for (int i = 0; i < n; ++i) {
          fgets(line, sizeof(line), fp);
          sscanf(line, "%d %d %d %d",
                 &z->xptidx2[i], &z->seplineidx2[i],
                 &z->segmidx2[i], &z->reverse_segm2[i]);
        }
      }
    }
  }

  fclose(fp);
  printf("Finish loading GridZoneInfo %s from file %s\n", z->name ? z->name : "Unnamed", filename);
  return z;
}
void print_GridZoneInfo(GridZoneInfo* gridzoneinfo)
{
  if (!gridzoneinfo) {
    printf("GridZoneInfo is NULL.\n");
    return;
  }

  printf("== GridZoneInfo Basic Info ==\n");
  printf("Topo: %s\n",  gridzoneinfo->topo ? gridzoneinfo->topo : "(null)");
  printf("Name: %s\n",  gridzoneinfo->name ? gridzoneinfo->name : "(null)");
  printf("nr:   %d\n",  gridzoneinfo->nr);

  printf("\n== Tracing Start Points ==\n");
  for (int i = 0; i < gridzoneinfo->nr; ++i) {
    printf("[%d] R=%.12f  Z=%.12f  guard_start=%.12f  guard_end=%.12f  pasmin=%.12f\n",
           i,
           gridzoneinfo->start_point_R[i],
           gridzoneinfo->start_point_Z[i],
           gridzoneinfo->guard_start[i],
           gridzoneinfo->guard_end[i],
           gridzoneinfo->pasmin[i]);
  }

  if (gridzoneinfo->start_curve && gridzoneinfo->start_curve->points) {
    printf("\n== Start Curve (%zu points) ==\n", gridzoneinfo->start_curve->n_point);
    for (size_t i = 0; i < gridzoneinfo->start_curve->n_point; ++i) {
      printf("[%zu] R=%.12f  Z=%.12f\n",
             i,
             gridzoneinfo->start_curve->points[i].x,
             gridzoneinfo->start_curve->points[i].y);
    }
  } else {
    printf("\n== Start Curve: NULL ==\n");
  }

  if (gridzoneinfo->end_curve && gridzoneinfo->end_curve->points) {
    printf("\n== End Curve (%zu points) ==\n", gridzoneinfo->end_curve->n_point);
    for (size_t i = 0; i < gridzoneinfo->end_curve->n_point; ++i) {
      printf("[%zu] R=%.12f  Z=%.12f\n",
             i,
             gridzoneinfo->end_curve->points[i].x,
             gridzoneinfo->end_curve->points[i].y);
    }
  } else {
    printf("\n== End Curve: NULL ==\n");
  }

  printf("\n== First Poloidal Segment Index Info (%d segments) ==\n", gridzoneinfo->n_polsegm1);
  for (int i = 0; i < gridzoneinfo->n_polsegm1; ++i) {
    printf("  Segment %d: xptidx1=%d seplineidx1=%d segmidx1=%d reverse=%d\n",
           i,
           gridzoneinfo->xptidx1[i],
           gridzoneinfo->seplineidx1[i],
           gridzoneinfo->segmidx1[i],
           gridzoneinfo->reverse_segm1[i]);
  }

  if (gridzoneinfo->n_polsegm2 != -1) {
    printf("\n== Second Poloidal Segment Index Info (%d segments) ==\n", gridzoneinfo->n_polsegm2);
    for (int i = 0; i < gridzoneinfo->n_polsegm2; ++i) {
      printf("  Segment %d: xptidx2=%d seplineidx2=%d segmidx2=%d reverse=%d\n",
             i,
             gridzoneinfo->xptidx2[i],
             gridzoneinfo->seplineidx2[i],
             gridzoneinfo->segmidx2[i],
             gridzoneinfo->reverse_segm2[i]);
    }
  } else {
    printf("\n== Second Poloidal Segment Index Info: Not Available (-1) ==\n");
  }
}