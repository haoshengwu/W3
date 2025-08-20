#include "curve.h"
#include <stdio.h>
#include <math.h>


#define TOL_CURVE 1.0E-10 // Tolerance for in-segment comparison
#define EPS_CURVE 1.0E-12 // Tolerance for in-segment comparison

OldCurve* create_oldcurve(size_t n_point)
{   
  if(!n_point)
  {
    fprintf(stderr, "Error: n_point is 0!\n");
    return NULL;
  }
  OldCurve* c = malloc(sizeof(OldCurve));
  if (!c)
  {
    fprintf(stderr, "Error: Memory allocation failed for a curve!\n");
    return NULL;
  } 
  c->n_point = n_point;
  c->points = allocate_2d_array(n_point,2);
  return c;
}

void free_oldcurve(OldCurve* c) {
    if (!c) return;
    free_2d_array(c->points);  // Free row pointers
    free(c);         // Free the OldCurve struct
}


int copy_oldcurve(OldCurve* c1, OldCurve* c2)
{
    if (!c1 || !c2) {
        fprintf(stderr, "Null pointer passed to copy_oldcurve.\n");
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

Curve* create_curve(size_t init_capacity) 
{
    if (init_capacity == 0) init_capacity = 128;

    Curve* c = malloc(sizeof(Curve));
    if (!c) 
    {
        fprintf(stderr, "Error: failed to allocate Curve structure.\n");
        return NULL;
    }

    c->points = malloc(init_capacity * sizeof(CurvePoint));
    if (!c->points) 
    {
        fprintf(stderr, "Error: failed to allocate points array.\n");
        free(c);
        return NULL;
    }

    c->n_point = 0;
    c->capacity = init_capacity;
    return c;
}

double get_curve_x(Curve* c, size_t idx)
{
  return c->points[idx].x;
}

double get_curve_y(Curve* c, size_t idx)
{
  return c->points[idx].y;
}

size_t get_curve_npt(Curve* c)
{
  return c->n_point;
}

void expand_curve_size_with_NaN(Curve* c, size_t size)
{
  if(!c||!size)
  {
    fprintf(stderr,"Empty inputs for expand_curve_size_with_NaN.\n");
    exit(EXIT_FAILURE);
  }

  for(int i=0;i<size;i++)
  {
    add_last_point_curve(c, NAN, NAN);
  }
}

void free_curve(Curve* c) 
{
    if (c) 
    {
        free(c->points);
        c->points=NULL;
        free(c);
    }
}

int add_last_point_curve(Curve* c, double x, double y) 
{
    if (!c) 
    {
      fprintf(stderr,"Empty input for add_point.\n");
      return -1;
    }
    // Resize array if capacity is full
    if (c->n_point == c->capacity) {
        size_t new_capacity = c->capacity * 2;
        CurvePoint* new_points = realloc(c->points, new_capacity * sizeof(CurvePoint));
        if (!new_points) 
        {
            fprintf(stderr, "Error: realloc failed while adding point.\n");
            return -1;
        }
        c->points = new_points;
        c->capacity = new_capacity;
    }
    c->points[c->n_point].x = x;
    c->points[c->n_point].y = y;
    c->n_point++;
    return 0;
}

int set_point_curve(Curve* c, size_t idx, double x, double y) 
{
    if (!c || idx >= c->n_point) 
    {
      fprintf(stderr,"Empty input for set_point_curve.\n");
      return -1;
    }
    c->points[idx].x = x;
    c->points[idx].y = y;
    return 0;
}

int delete_last_point(Curve* c) {
    if (!c) 
    {
        fprintf(stderr, "Error: NULL curve passed to delete_last_point.\n");
        return -1;
    }
    if (c->n_point == 0) 
    {
        fprintf(stderr, "Warning: delete_last_point called on empty curve.\n");
        return -1;
    }

    c->n_point--;  // Just decrease logical length
    return 0;
}

void write_curve(const char *filename, const Curve *c) 
{
    if (!filename || !c) {
        fprintf(stderr, "Error: NULL input to write_curve.\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot open file \"%s\" \n",filename);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < c->n_point; ++i) {
        if (fprintf(fp, "%.12f %.12f\n", c->points[i].x, c->points[i].y) < 0) 
        {
            fprintf(stderr, "Error: write failed while writing to \"%s\"\n", filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
    }
    printf("Finish writing curev to %s file.\n",filename);
    fclose(fp);
}

void print_curve(const Curve *c)
{
    if (!c) 
    {
        fprintf(stderr, "Error: NULL pointer passed to print_curve.\n");
        return;
    }

    for (size_t i = 0; i < c->n_point; ++i) 
    {
        printf("%zu: %.15g %.15g\n",
               i, c->points[i].x, c->points[i].y);
    }
}

Curve* copy_curve(Curve* c)
{
  // Check for null input
  if (!c) {
    fprintf(stderr, "Error: NULL input to copy_curve.\n");
    exit(EXIT_FAILURE);
  }

  // Defensive check: ensure number of points does not exceed capacity
  if (c->n_point > c->capacity) {
    fprintf(stderr, "Error: Curve n_point (%zu) exceeds capacity (%zu). Possible memory corruption.\n",
            c->n_point, c->capacity);
    exit(EXIT_FAILURE);
  }

  // If the curve has points, its points array must not be NULL
  if (c->n_point > 0 && !c->points) {
    fprintf(stderr, "Error: Curve has n_point > 0 but points array is NULL.\n");
    exit(EXIT_FAILURE);
  }

  // Allocate a new curve with the same capacity
  Curve* new_c = create_curve(c->capacity);
  if (!new_c) {
    fprintf(stderr, "Error: Failed to allocate new Curve in copy_curve.\n");
    exit(EXIT_FAILURE);
  }

  // Copy used points only (not the entire capacity)
  if (c->n_point > 0) {
    memcpy(new_c->points, c->points, c->n_point * sizeof(CurvePoint));
    new_c->n_point = c->n_point;
  }
  return new_c;
}


Curve* convert_ddl_to_curve(DLListNode* head)
{
  if(!head)
  {
    fprintf(stderr, "Empty inputs for convert_ddl_to_curve.\n");
    exit(EXIT_FAILURE);
  }
  Curve* curve = create_curve(0);
  DLListNode* cur_node = head;
  while(cur_node)
  {
    if(add_last_point_curve(curve, cur_node->r, cur_node->z)!=0)
    {
      fprintf(stderr, "Unexpected error: cannot add point to the curve in convert_ddl_to_curve\n");
      exit(EXIT_FAILURE);
    }
    cur_node=cur_node->next;
  }
  if(curve->n_point<2)
  {
    fprintf(stderr, "The curve only has %zu point(s) in convert_ddl_to_curve\n", curve->n_point);
    exit(EXIT_FAILURE);
  }
  return curve;
}

/* 2D cross product: (B - A) × (C - A) */
static inline double cross(double ax, double ay, 
                           double bx, double by,
                           double cx, double cy) 
{
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

/* Check if (cx, cy) lies on segment (ax, ay)-(bx, by) */
static inline int on_segment(double ax, double ay, 
                             double bx, double by, 
                             double cx, double cy) 
{
    return fmin(ax, bx) - EPS_CURVE <= cx && cx <= fmax(ax, bx) + EPS_CURVE &&
           fmin(ay, by) - EPS_CURVE <= cy && cy <= fmax(ay, by) + EPS_CURVE;
}

/* Return 0 if segments intersect (including endpoints), 1 if not */
int has_intersection(double x1, double y1, double x2, double y2,
                     double x3, double y3, double x4, double y4)
{
    // Step 1: Bounding box check
    if (fmax(x1, x2) + EPS_CURVE < fmin(x3, x4) ||
        fmax(x3, x4) + EPS_CURVE < fmin(x1, x2) ||
        fmax(y1, y2) + EPS_CURVE < fmin(y3, y4) ||
        fmax(y3, y4) + EPS_CURVE < fmin(y1, y2)) 
    {
        return 1;
    }

    // Step 2: Cross product signs
    double d1 = cross(x3, y3, x4, y4, x1, y1);
    double d2 = cross(x3, y3, x4, y4, x2, y2);
    double d3 = cross(x1, y1, x2, y2, x3, y3);
    double d4 = cross(x1, y1, x2, y2, x4, y4);

    if ((d1 * d2 < 0) && (d3 * d4 < 0)) return 0;

    if (fabs(d1) < EPS_CURVE && on_segment(x3, y3, x4, y4, x1, y1)) return 0;
    if (fabs(d2) < EPS_CURVE && on_segment(x3, y3, x4, y4, x2, y2)) return 0;
    if (fabs(d3) < EPS_CURVE && on_segment(x1, y1, x2, y2, x3, y3)) return 0;
    if (fabs(d4) < EPS_CURVE && on_segment(x1, y1, x2, y2, x4, y4)) return 0;

    return 1;
}

int get_intersection_point(double x1, double y1, double x2, double y2,
                           double x3, double y3, double x4, double y4,
                           double* intsect_x, double* intsect_y)
{
    // First check if they intersect
    // Bounding box rejection
    if (fmax(x1, x2) + EPS_CURVE < fmin(x3, x4) ||
        fmax(x3, x4) + EPS_CURVE < fmin(x1, x2) ||
        fmax(y1, y2) + EPS_CURVE < fmin(y3, y4) ||
        fmax(y3, y4) + EPS_CURVE < fmin(y1, y2)) 
    {
        *intsect_x=NAN;
        *intsect_y=NAN;
        return 1;
    }

    // Cross products
    double d1 = cross(x3, y3, x4, y4, x1, y1);
    double d2 = cross(x3, y3, x4, y4, x2, y2);
    double d3 = cross(x1, y1, x2, y2, x3, y3);
    double d4 = cross(x1, y1, x2, y2, x4, y4);

    // General case
    if ((d1 * d2 < 0) && (d3 * d4 < 0)) {
        // Use parametric intersection: compute t
        double dx1 = x2 - x1;
        double dy1 = y2 - y1;
        double dx2 = x4 - x3;
        double dy2 = y4 - y3;

        double denom = dx1 * dy2 - dy1 * dx2;
        if (fabs(denom) < EPS_CURVE) 
        {
            *intsect_x=NAN;
            *intsect_y=NAN;
            return 1; // Parallel or degenerate
        }
        double t = ((x3 - x1) * dy2 - (y3 - y1) * dx2) / denom;

        *intsect_x = x1 + t * dx1;
        *intsect_y = y1 + t * dy1;
        return 0;
    }

    // Special cases: endpoints overlap
    if (fabs(d1) < EPS_CURVE && fmin(x3,x4)-EPS_CURVE <= x1 && x1 <= fmax(x3,x4)+EPS_CURVE &&
                                 fmin(y3,y4)-EPS_CURVE <= y1 && y1 <= fmax(y3,y4)+EPS_CURVE) 
    {
        *intsect_x = x1; *intsect_y = y1; return 0;
    }
    if (fabs(d2) < EPS_CURVE && fmin(x3,x4)-EPS_CURVE <= x2 && x2 <= fmax(x3,x4)+EPS_CURVE &&
                                 fmin(y3,y4)-EPS_CURVE <= y2 && y2 <= fmax(y3,y4)+EPS_CURVE) 
    {
        *intsect_x = x2; *intsect_y = y2; return 0;
    }
    if (fabs(d3) < EPS_CURVE && fmin(x1,x2)-EPS_CURVE <= x3 && x3 <= fmax(x1,x2)+EPS_CURVE &&
                                 fmin(y1,y2)-EPS_CURVE <= y3 && y3 <= fmax(y1,y2)+EPS_CURVE) 
    {
        *intsect_x = x3; *intsect_y = y3; return 0;
    }
    if (fabs(d4) < EPS_CURVE && fmin(x1,x2)-EPS_CURVE <= x4 && x4 <= fmax(x1,x2)+EPS_CURVE &&
                                 fmin(y1,y2)-EPS_CURVE <= y4 && y4 <= fmax(y1,y2)+EPS_CURVE) 
    {
        *intsect_x = x4; *intsect_y = y4; return 0;
    }
    return 1;
}
/*************************************************
 * Use for TwoDim Grid Generation
**************************************************/
double total_length_curve(const Curve *c)
{
    if (!c || c->n_point < 2) return 0.0;

    double length = 0.0;

    for (size_t i = 1; i < c->n_point; ++i) {
        double dx = c->points[i].x - c->points[i - 1].x;
        double dy = c->points[i].y - c->points[i - 1].y;
        /* hypot(dx, dy) = sqrt(dx*dx + dy*dy) with better numerical stability */
        length += hypot(dx, dy);
    }
    return length;
}

double length_curve(Curve *curve, size_t ith)
{
    if (!curve) {
        fprintf(stderr, "Error: NULL input to length_curve.\n");
        exit(EXIT_FAILURE);
    }
    if (ith > curve->n_point) {
        fprintf(stderr, "Error: n (%zu) is larger than curve->n_point (%zu).\n", ith, curve->n_point);
        exit(EXIT_FAILURE);
    }
    if (ith < 2) {
        return 0.0;  // 0 or 1 point → no length
    }

    double length = 0.0;
    for (size_t i = 1; i < ith; ++i) {
        double dx = curve->points[i].x - curve->points[i - 1].x;
        double dy = curve->points[i].y - curve->points[i - 1].y;
        length += hypot(dx, dy);
    }
    return length;
}

int indcrb_curve(Curve* curve, CurvePoint* point, double d)
{
    if (!curve || !point || curve->n_point < 2) {
        fprintf(stderr, "Error: invalid input to indcrb_CARRE.\n");
        exit(EXIT_FAILURE);
    }

    double mumin = 1.0;
    int indcrb = 0;
    double ax, ay, bx, by, dist, mu;
    double dd = 0.0;

    // Step 1: check if point is very close to the first point
    ax = curve->points[0].x - point->x;
    ay = curve->points[0].y - point->y;
    dist = sqrt(ax * ax + ay * ay);

    double dx1 = curve->points[1].x - curve->points[0].x;
    double dy1 = curve->points[1].y - curve->points[0].y;
    double first_segment_length = hypot(dx1, dy1);

    if (dist < TOL_CURVE && first_segment_length > d) {
        return 0;
    }

    for (size_t i = 0; i < curve->n_point - 1; ++i) {
        double dx = curve->points[i + 1].x - curve->points[i].x;
        double dy = curve->points[i + 1].y - curve->points[i].y;
        double seg_len = hypot(dx, dy);
        dd += seg_len;

        if (dd >= d) {
            ax = curve->points[i].x - point->x;
            ay = curve->points[i].y - point->y;
            bx = curve->points[i + 1].x - point->x;
            by = curve->points[i + 1].y - point->y;

            dist = hypot(bx, by);
            if (dist < TOL_CURVE) 
            {
                return (i < curve->n_point - 2) ? (int)(i + 1) : (int)i;
            } else 
            {
                double norm_a = hypot(ax, ay);
                double norm_b = hypot(bx, by);
                mu = (ax * bx + ay * by) / (norm_a * norm_b);
                if (mu < mumin) 
                 {
                    mumin = mu;
                    indcrb = (int)i;
                }
            }
        }
    }
    return indcrb;
}

double ruban_curve(Curve* curve, CurvePoint* point, double d)
{
    if (!curve || !point || curve->n_point < 2) {
        fprintf(stderr, "Invalid input to ruban_CARRE_Curve.\n");
        exit(EXIT_FAILURE);
    }

    double ruban = 0.0; 
    double dist = 0.0;  
    int ind;
    double dx, dy;

    // Find which segment the point lies on (using d to help guide the search)
    ind = indcrb_curve(curve, point, d);

    dist = length_curve(curve, ind);

    dx = point->x - curve->points[ind].x;
    dy = point->y - curve->points[ind].y;
    ruban = dist + hypot(dx, dy);

    return ruban;
}

void coordnates_in_curve(Curve *curve, double d, CurvePoint *point) {
    if (!curve || curve->n_point < 2 || !point) {
        fprintf(stderr, "Error: invalid input to coord_in_curve.\n");
        exit(EXIT_FAILURE);
    }

    // Step 1: compute total curve length
    double total_len = 0.0;
    for (size_t i = 1; i < curve->n_point; ++i) {
        double dx = curve->points[i].x - curve->points[i - 1].x;
        double dy = curve->points[i].y - curve->points[i - 1].y;
        total_len += hypot(dx, dy);
    }

    // Step 2: handle endpoints with TOL_CURVEerance
    if (d <= TOL_CURVE) {
        point->x = curve->points[0].x;
        point->y = curve->points[0].y;
        return;
    }
    if (fabs(d - total_len) <= TOL_CURVE) {
        point->x = curve->points[curve->n_point - 1].x;
        point->y = curve->points[curve->n_point - 1].y;
        return;
    }

    if (d < 0.0 || d > total_len + TOL_CURVE) {
        fprintf(stderr, "Error: distance %.15g out of bounds [0, %.15g]\n", d, total_len);
        exit(EXIT_FAILURE);
    }

    // Step 3: walk through segments to locate d
    double accum = 0.0;

    for (size_t i = 0; i < curve->n_point - 1; ++i) {
        double x0 = curve->points[i].x;
        double y0 = curve->points[i].y;
        double x1 = curve->points[i + 1].x;
        double y1 = curve->points[i + 1].y;

        double dx = x1 - x0;
        double dy = y1 - y0;
        double seg_len = hypot(dx, dy);
        double next_accum = accum + seg_len;

        // Floating-point safe inclusion test
        if (d >= accum - TOL_CURVE && d <= next_accum + TOL_CURVE) {
            double remain = d - accum;
            if (remain <= TOL_CURVE) {
                point->x = curve->points[i].x;  // Near segment start
                point->y = curve->points[i].y;
            } else if (fabs(remain - seg_len) <= TOL_CURVE) {
                point->x = curve->points[i + 1].x;  // Near segment end
                point->x = curve->points[i + 1].y;  // Near segment end
            } else {
                double ratio = remain / seg_len;
                point->x = x0 + ratio * dx;
                point->y = y0 + ratio * dy;
            }
            return;
        }
        accum = next_accum;
    }

    // Should never reach here
    fprintf(stderr, "coord_in_curve: interpolation failure at d = %.15g\n", d);
    exit(EXIT_FAILURE);
}

//create a new curve which is due the connection of max to five curve
//REVERSE FIRST, THEN SKIP
Curve* connect_curves_for_curve(CurveWithOptions* list, int n_curve)
{
  if (!list || n_curve <= 0) return NULL;

  Curve* result = create_curve(256);
  if (!result) return NULL;

  for (int i = 0; i < n_curve; i++) {
    Curve* c = list[i].curve;
    if (!c || c->n_point == 0) {
      fprintf(stderr, "Error: segment %d is empty or null.\n", i);
      exit(EXIT_FAILURE);
    }

    // Determine iteration direction and index range
    int step  = list[i].reverse ? -1 : 1;
    int first = list[i].reverse ? c->n_point - 1 : 0;
    int last  = list[i].reverse ? 0 : c->n_point - 1;

    // Validate index range
    if ((step == 1 && first > last) ||
        (step == -1 && first < last)) {
      fprintf(stderr, "Error: segment %d has invalid index range after reverse.\n", i);
      exit(EXIT_FAILURE);
    }

    // Automatically skip first point if it overlaps with previous segment's last point
    int adjusted_first = first;
    if (result->n_point > 0) {
      CurvePoint* prev = &result->points[result->n_point - 1];
      CurvePoint* curr_first = &c->points[first];
      if (fabs(prev->x - curr_first->x) < 1e-10 &&
          fabs(prev->y - curr_first->y) < 1e-10) {
        adjusted_first += step;
      }
    }

    // After skipping overlap, check if segment is empty
    if ((step == 1 && adjusted_first > last) ||
        (step == -1 && adjusted_first < last)) {
      fprintf(stderr, "Error: segment %d has no valid points after automatic overlap skip.\n", i);
      exit(EXIT_FAILURE);
    }

    // Append valid points to the result curve
    for (int idx = adjusted_first; idx != last + step; idx += step) {
      add_last_point_curve(result, c->points[idx].x, c->points[idx].y);
    }
  }
  return result;
}