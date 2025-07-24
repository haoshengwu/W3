/**
 * @file    config.h
 * @brief   Configuration parameters for the project
 * @date    2024-01-10
 * 
 * This header contains all compile-time configuration constants.
 * Some values can be overridden during compilation.
 */

#ifndef CONFIG_H
#define CONFIG_H

/* ========== Build Configuration ========== */


/* ========== Numerical Constants ========== */
#define MAX_NUM_TRACING 40000 
#define NRELAX 3000

#define PI 3.14159265358979323846


/* ========== Precision Settings ========== */
#define EPSILON 1.0E-12

/* ========== 2DGRID and 3DGRID Settings ========== */
#define DEFAULT_POL_MARGIN 20
#define DEFAULT_RAD_MARGIN 20
#define DEFAULT_TOR_MARGIN 0


//for the inner and outer lege, the original resolution minus nlast/nfirst should larger than MIN_RESOLTUTION
#define MIN_RESOLTUTION 5

#define GRID_ALIGNMENT 32


#endif /* CONFIG_H */