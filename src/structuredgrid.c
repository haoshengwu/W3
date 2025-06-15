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
  z->npoint=-1;
  z->nr=-1;
  
  z->grid_curveset=NULL;

  z->start_point_R=NULL;
  z->start_point_Z=NULL;
  z->guard_start=NULL;
  z->guard_end=NULL;
  z->pasmin = NULL;
  z->norm_pol_dist=NULL;
 
  z->first_boundary=NULL;
//  z->second_boundary=NULL;

  z->end_curve=NULL;

  return z;
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
  // if((*z)->second_boundary)
  // {
  //   free_curve((*z)->second_boundary);
  //   (*z)->second_boundary=NULL;
  // }

  free_curve((*z)->end_curve);
  (*z)->end_curve=NULL;

  free(*z); 
  *z = NULL;
}

void write_input_from_GridZone(GridZone* gridzone, const char* filename)
{
    if (!gridzone || !filename) {
        fprintf(stderr, "Invalid input to write_input_from_GridZone.\n");
        exit(EXIT_FAILURE);
    }

    FILE* fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Invalid input to write_input_from_GridZone.\n");
        exit(EXIT_FAILURE);
    }

    // --- 1. Basic information ---
    fprintf(fp, "#name\n%s\n", gridzone->name ? gridzone->name : "Unnamed");
    fprintf(fp, "#nr\n%d\n", gridzone->nr);

    // --- 2. Start & end tracing info ---
    fprintf(fp, "#start_point_R Start_point_Z guardstart guardend pasmin\n");
    for (int i = 0; i < gridzone->nr; ++i) {
        fprintf(fp, "%.10f %.10f %.10f %.10f %.10f\n",
                gridzone->start_point_R[i],
                gridzone->start_point_Z[i],
                gridzone->guard_start[i],
                gridzone->guard_end[i],
                gridzone->pasmin[i]);
    }

    // --- first line distribution ---
    fprintf(fp, "#npoint\n%d\n", gridzone->npoint);
    fprintf(fp, "#norm_pol_dist\n");
    for (int i = 0; i < gridzone->npoint; ++i) {
        fprintf(fp, "%.10f\n", gridzone->norm_pol_dist[i]);
    }

    // --- 3. Boundary info ---
    fprintf(fp, "#first_boundary\n");
    if (gridzone->first_boundary && gridzone->first_boundary->points) {
        fprintf(fp, "%zu\n", gridzone->first_boundary->n_point);  // Number of points
        for (int i = 0; i < gridzone->first_boundary->n_point; ++i) {
            fprintf(fp, "%.10f %.10f\n",
                    gridzone->first_boundary->points[i][0],
                    gridzone->first_boundary->points[i][1]);
        }
    } else {
        fprintf(fp, "0\n# No first_boundary defined.\n");
    }

    // --- 4. Target info ---
    fprintf(fp, "#end_curve\n");
    if (gridzone->end_curve && gridzone->end_curve->points) {
        fprintf(fp, "%zu\n", gridzone->end_curve->n_point);  // Number of points
        for (int i = 0; i < gridzone->end_curve->n_point; ++i) {
            fprintf(fp, "%.10f %.10f\n",
                    gridzone->end_curve->points[i][0],
                    gridzone->end_curve->points[i][1]);
        }
    } else {
        fprintf(fp, "0\n# No end_curve defined.\n");
    }
    fclose(fp);
    printf("Finish writing the %s to file %s\n",gridzone->name, filename);
}

GridZone* load_GridZone_from_input(const char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open input file");
        exit(EXIT_FAILURE);
    }

    GridZone* z = allocate_GridZone();

    char line[256];

    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "#name", 5) == 0) {
            fgets(line, sizeof(line), fp);
            line[strcspn(line, "\r\n")] = 0;
            size_t len = strlen(line);
            z->name = malloc(len + 1);
            if (!z->name) {
                fprintf(stderr, "Memory allocation failed for name.\n");
                exit(EXIT_FAILURE);
            }
            strcpy(z->name, line);
        }

        else if (strncmp(line, "#nr", 3) == 0) {
            fgets(line, sizeof(line), fp);
            sscanf(line, "%d", &z->nr);
            int nr = z->nr;
            z->start_point_R = malloc(nr * sizeof(double));
            z->start_point_Z = malloc(nr * sizeof(double));
            z->guard_start   = malloc(nr * sizeof(double));
            z->guard_end     = malloc(nr * sizeof(double));
            z->pasmin        = malloc(nr * sizeof(double));
            if (!z->start_point_R || !z->start_point_Z || !z->guard_start ||
                !z->guard_end || !z->pasmin) {
                fprintf(stderr, "Memory allocation failed for tracing arrays.\n");
                exit(EXIT_FAILURE);
            }
        }

        else if (strncmp(line, "#start_point_R", 14) == 0) {
            for (int i = 0; i < z->nr; ++i) {
                fgets(line, sizeof(line), fp);
                sscanf(line, "%lf %lf %lf %lf %lf",
                       &z->start_point_R[i],
                       &z->start_point_Z[i],
                       &z->guard_start[i],
                       &z->guard_end[i],
                       &z->pasmin[i]);
            }
        }

        else if (strncmp(line, "#npoint", 7) == 0) {
            fgets(line, sizeof(line), fp);
            sscanf(line, "%d", &z->npoint);
            z->norm_pol_dist = malloc(z->npoint * sizeof(double));
            if (!z->norm_pol_dist) {
                fprintf(stderr, "Memory allocation failed for norm_pol_dist.\n");
                exit(EXIT_FAILURE);
            }
        }

        else if (strncmp(line, "#norm_pol_dist", 14) == 0) {
            for (int i = 0; i < z->npoint; ++i) {
                fgets(line, sizeof(line), fp);
                sscanf(line, "%lf", &z->norm_pol_dist[i]);
            }
        }

        else if (strncmp(line, "#first_boundary", 15) == 0) {
            fgets(line, sizeof(line), fp);
            size_t n;
            sscanf(line, "%zu", &n);
            if (n > 0) {
                z->first_boundary = create_curve(n);
                if (!z->first_boundary) {
                    fprintf(stderr, "Memory allocation failed for first_boundary.\n");
                    exit(EXIT_FAILURE);
                }
                for (size_t i = 0; i < n; ++i) {
                    fgets(line, sizeof(line), fp);
                    sscanf(line, "%lf %lf",
                           &z->first_boundary->points[i][0],
                           &z->first_boundary->points[i][1]);
                }
            }
        }

        else if (strncmp(line, "#end_curve", 10) == 0) {
            fgets(line, sizeof(line), fp);
            size_t n;
            sscanf(line, "%zu", &n);
            if (n > 0) {
                z->end_curve = create_curve(n);
                if (!z->end_curve) {
                    fprintf(stderr, "Memory allocation failed for end_curve.\n");
                    exit(EXIT_FAILURE);
                }
                for (size_t i = 0; i < n; ++i) {
                    fgets(line, sizeof(line), fp);
                    sscanf(line, "%lf %lf",
                           &z->end_curve->points[i][0],
                           &z->end_curve->points[i][1]);
                }
            }
        }
    }
    fclose(fp);
    printf("Finish loading %s from file %s\n",z->name, filename);
    return z;
}

void print_GridZone(GridZone* gridzone)
{
    if (!gridzone) {
        printf("GridZone is NULL.\n");
        return;
    }

    // --- 1. Basic information ---
    printf("== GridZone Basic Info ==\n");
    printf("Name: %s\n", gridzone->name ? gridzone->name : "(null)");
    printf("nr: %d\n", gridzone->nr);
    printf("npoint: %d\n", gridzone->npoint);  // ✅ 显式打印 npoint

    // --- 2. Start & end tracing info ---
    printf("\n== Tracing Start Points ==\n");
    for (int i = 0; i < gridzone->nr; ++i) {
        printf("[%d] R=%.10f  Z=%.10f  guard_start=%.10f  guard_end=%.10f  pasmin=%.10f\n",
               i,
               gridzone->start_point_R[i],
               gridzone->start_point_Z[i],
               gridzone->guard_start[i],
               gridzone->guard_end[i],
               gridzone->pasmin[i]);
    }

    // --- 3. norm_pol_dist ---
    printf("\n== Norm Poloidal Distribution (%d points) ==\n", gridzone->npoint);
    for (int i = 0; i < gridzone->npoint; ++i) {
        printf("[%d] %.10f\n", i, gridzone->norm_pol_dist[i]);
    }

    // --- 4. End curve ---
    if (gridzone->end_curve && gridzone->end_curve->points) {
        printf("\n== End Curve (%zu points) ==\n", gridzone->end_curve->n_point);
        for (size_t i = 0; i < gridzone->end_curve->n_point; ++i) {
            printf("[%zu] R=%.10f  Z=%.10f\n",
                   i,
                   gridzone->end_curve->points[i][0],
                   gridzone->end_curve->points[i][1]);
        }
    } else {
        printf("\n== End Curve: NULL ==\n");
    }
}
