#include "polsegminfo.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

PolSegmsInfo* create_PolSegmsInfo(int n_polsegms)
{
  if (n_polsegms <= 0) return NULL;

  PolSegmsInfo* info = (PolSegmsInfo*)malloc(sizeof(PolSegmsInfo));
  if (!info) {
    fprintf(stderr, "Failed to allocate PolSegmsInfo.\n");
    exit(EXIT_FAILURE);
  }

  info->topo = NULL;
  info->n_polsegms = n_polsegms;
  info->polsegments = (PolSegStr**)malloc(n_polsegms * sizeof(PolSegStr*));
  if (!info->polsegments) {
    fprintf(stderr, "Failed to allocate PolSegStr* array.\n");
    free(info);
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < n_polsegms; ++i)
    info->polsegments[i] = NULL;

  info->xptidx=malloc(n_polsegms*sizeof(int));
  info->seplineidx=malloc(n_polsegms*sizeof(int));
  info->segmidx=malloc(n_polsegms*sizeof(int));
  info->reverse_segm=malloc(n_polsegms*sizeof(int));

  if(!info->xptidx||!info->seplineidx
    ||!info->segmidx||!info->reverse_segm)
  {
    fprintf(stderr, "Failed to allocate Entries of PolSegmsInfo.\n");
    exit(EXIT_FAILURE);
  }
  return info;
}

PolSegStr* create_PolSegStr(int n_points, char* name)
{
  if (n_points <= 0 || !name) return NULL;

  PolSegStr* seg = (PolSegStr*)malloc(sizeof(PolSegStr));
  if (!seg) {
    fprintf(stderr, "Failed to allocate PolSegStr.\n");
    exit(EXIT_FAILURE);
  }

  seg->n_points = n_points;
  seg->norm_dist = (double*)malloc(n_points * sizeof(double));
  if (!seg->norm_dist) {
    fprintf(stderr, "Failed to allocate norm_dist.\n");
    free(seg);
    exit(EXIT_FAILURE);
  }

  size_t len = strcspn(name, "\r\n");
  seg->name = (char*)malloc(len + 1);
  if (!seg->name) {
    fprintf(stderr, "Failed to allocate name.\n");
    free(seg->norm_dist);
    free(seg);
    exit(EXIT_FAILURE);
  }
  strncpy(seg->name, name, len);
  seg->name[len] = '\0';

  for (int i = 0; i < n_points; ++i)
    seg->norm_dist[i] = NAN;

  return seg;
}

void free_PolSegStr(PolSegStr* polseg)
{
  if (!polseg) return;
  free(polseg->name);
  free(polseg->norm_dist);
  free(polseg);
}

void free_PolSegmsInfo(PolSegmsInfo* polseginfo)
{
  if (!polseginfo) return;

  free(polseginfo->topo);

  if (polseginfo->polsegments) 
  {
    for (int i = 0; i < polseginfo->n_polsegms; ++i) 
    {
      if (polseginfo->polsegments[i]) 
      {
        free_PolSegStr(polseginfo->polsegments[i]);
        polseginfo->polsegments[i] = NULL;
      }
    }
    free(polseginfo->polsegments);
    polseginfo->polsegments = NULL;
  }

  free(polseginfo->xptidx);
  polseginfo->xptidx=NULL;
  free(polseginfo->seplineidx);
  polseginfo->seplineidx=NULL;
  free(polseginfo->segmidx);
  polseginfo->segmidx=NULL;
  free(polseginfo->reverse_segm);
  polseginfo->reverse_segm=NULL;
  free(polseginfo);
}

void write_PolSegStr(PolSegStr* polseg, FILE* fp)
{
  if (!polseg || !fp) {
    fprintf(stderr, "Invalid input to write_PolSegStr.\n");
    return;
  }

  fprintf(fp, "#segment\n%s\n", polseg->name ? polseg->name : "UnnamedSegment");
  fprintf(fp, "#size\n%d\n", polseg->n_points);
  fprintf(fp, "#normal distribution\n");
  for (int i = 0; i < polseg->n_points; ++i) {
    fprintf(fp, "%.12f\n", polseg->norm_dist[i]);
  }
}

void write_PolSegmsInfo(PolSegmsInfo* polseginfo, const char* filename)
{
  if (!polseginfo || !filename) {
    fprintf(stderr, "Invalid input to write_PolSegmsInfo.\n");
    return;
  }

  FILE* fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s for writing.\n", filename);
    exit(EXIT_FAILURE);
  }

  fprintf(fp, "#topo\n%s\n", polseginfo->topo ? polseginfo->topo : "NaN");
  fprintf(fp, "#nsegments\n%d\n", polseginfo->n_polsegms);

  for (int i = 0; i < polseginfo->n_polsegms; ++i) 
  {
    if (polseginfo->polsegments[i]) {
      write_PolSegStr(polseginfo->polsegments[i], fp);
    }
  }

  fprintf(fp, "#xptidx seplineidx segmidx reverse\n");
  for (int i = 0; i < polseginfo->n_polsegms; ++i) 
  {
    fprintf(fp, "%d %d %d %d\n", polseginfo->xptidx[i],
                                 polseginfo->seplineidx[i],
                                 polseginfo->segmidx[i], 
                                 polseginfo->reverse_segm[i]);
  }

  fclose(fp);
  printf("Finished writing PolSegmsInfo to file %s\n", filename);
}

void print_PolSegStr(PolSegStr* polseg)
{
  if (!polseg) {
    printf("  [Segment] (null)\n");
    return;
  }

  printf("  Segment: %s\n", polseg->name ? polseg->name : "(Unnamed)");
  printf("    Number of points: %d\n", polseg->n_points);
  printf("    Normalized distribution:\n");

  for (int i = 0; i < polseg->n_points; ++i) {
    printf("      [%3d] %.12f\n", i, polseg->norm_dist[i]);
  }
}

void print_PolSegmsInfo(PolSegmsInfo* polseginfo)
{
  if (!polseginfo) {
    printf("[PolSegmsInfo] (null)\n");
    return;
  }

  printf("=== Poloidal Segment Info ===\n");
  printf("Topo: %s\n", polseginfo->topo ? polseginfo->topo : "NaN");
  printf("Number of segments: %d\n", polseginfo->n_polsegms);

  for (int i = 0; i < polseginfo->n_polsegms; ++i) {
    printf("  [Segment %d]\n", i);
    print_PolSegStr(polseginfo->polsegments[i]);
  }
  
  printf("=== xptidx seplineidx segmidx reverse ===\n");
  for (int i = 0; i < polseginfo->n_polsegms; ++i) 
  {
    printf("  [Segment %d]  ", i);
    printf("%d %d %d %d\n", polseginfo->xptidx[i],
                                 polseginfo->seplineidx[i],
                                 polseginfo->segmidx[i], 
                                 polseginfo->reverse_segm[i]);
  }

  printf("=== End of Segment Info ===\n");
}

PolSegmsInfo* read_PolSegmsInfo_from_file(const char* filename)
{
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    perror("Failed to open input file");
    exit(EXIT_FAILURE);
  }

  char line[256];
  char* topo_tmp = NULL;
  int nsegments = -1;

  PolSegmsInfo* info = NULL;
  int current = -1;

  while (fgets(line, sizeof(line), fp)) {
    size_t len = strcspn(line, "\r\n");
    line[len] = '\0';

    if (strncmp(line, "#topo", 5) == 0) {
      if (!fgets(line, sizeof(line), fp)) goto error;
      len = strcspn(line, "\r\n");
      line[len] = '\0';

      topo_tmp = (char*)malloc(len + 1);
      if (!topo_tmp) {
        fprintf(stderr, "Failed to allocate memory for topo_tmp.\n");
        goto error;
      }
      memcpy(topo_tmp, line, len);
      topo_tmp[len] = '\0';

    } else if (strncmp(line, "#nsegments", 10) == 0) {
      if (!fgets(line, sizeof(line), fp)) goto error;
      len = strcspn(line, "\r\n");
      line[len] = '\0';

      if (sscanf(line, "%d", &nsegments) != 1 || nsegments <= 0) {
        fprintf(stderr, "Invalid segment count.\n");
        goto error;
      }

      info = create_PolSegmsInfo(nsegments);

      // copy topo to info->topo
      if (topo_tmp) {
        size_t topo_len = strlen(topo_tmp);
        info->topo = (char*)malloc(topo_len + 1);
        if (!info->topo) {
          fprintf(stderr, "Failed to allocate memory for info->topo.\n");
          goto error;
        }
        memcpy(info->topo, topo_tmp, topo_len);
        info->topo[topo_len] = '\0';
        free(topo_tmp);
        topo_tmp = NULL;
      }

    } else if (strncmp(line, "#segment", 8) == 0) {
      if (!info) goto error;

      // read segment name
      if (!fgets(line, sizeof(line), fp)) goto error;
      len = strcspn(line, "\r\n");
      line[len] = '\0';
      char segname[128];
      strncpy(segname, line, sizeof(segname) - 1);
      segname[sizeof(segname) - 1] = '\0';

      // read #size
      if (!fgets(line, sizeof(line), fp)) goto error;
      len = strcspn(line, "\r\n");
      line[len] = '\0';
      if (strncmp(line, "#size", 5) != 0) goto error;

      // read points
      if (!fgets(line, sizeof(line), fp)) goto error;
      len = strcspn(line, "\r\n");
      line[len] = '\0';
      int npts = 0;
      if (sscanf(line, "%d", &npts) != 1 || npts <= 0) goto error;

      // read #normal distribution
      if (!fgets(line, sizeof(line), fp)) goto error;
      len = strcspn(line, "\r\n");
      line[len] = '\0';
      if (strncmp(line, "#normal distribution", 20) != 0) goto error;

      ++current;
      if (current >= info->n_polsegms) goto error;

      PolSegStr* seg = create_PolSegStr(npts, segname);
      if (!seg) goto error;

      info->polsegments[current] = seg;

      // read distribution
      for (int i = 0; i < npts; ++i) {
        if (!fgets(line, sizeof(line), fp)) goto error;
        len = strcspn(line, "\r\n");
        line[len] = '\0';
        double val = 0.0;
        if (sscanf(line, "%lf", &val) != 1) goto error;
        seg->norm_dist[i] = val;
      }
    } 
    else if (strncmp(line, "#xptidx seplineidx segmidx reverse", 34) == 0)
    {
      for (int i = 0; i < info->n_polsegms; ++i) {
          if (!fgets(line, sizeof(line), fp)) goto error;
          if (sscanf(line, "%d %d %d %d", 
                     &info->xptidx[i], &info->seplineidx[i], 
                     &info->segmidx[i], &info->reverse_segm[i]) != 4) 
            {
              goto error;
            }
      }
    }
  }
  if (current + 1 != info->n_polsegms)
  {
    fprintf(stderr, "Warning: number of segments read (%d) does not match expected (%d).\n",
            current + 1, info->n_polsegms);
    goto error;
  }

  fclose(fp);
  return info;

error:
  fclose(fp);
  if (topo_tmp) free(topo_tmp);
  if (info) free_PolSegmsInfo(info);
  fprintf(stderr, "Error occurred while parsing file: %s\n", filename);
  exit(EXIT_FAILURE);
}