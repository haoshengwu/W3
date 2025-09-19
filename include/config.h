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
#define PI 3.14159265358979323846


/* ========== 2D GRID GENERATION ========== */

//FROM DIVGEO&CARRE:
//
#define TGT_GUARD_HEAD 0.16
#define TGT_GUARD_TAIL 0.16

#define PASMIN_CORE 1.0E-3
// #define PASMIN_SOL 1.0E-3 //High resolution
#define PASMIN_SOL 4.0E-3 //Low resolution


//#define PASMIN_SOL 2.0E-3 //High resolution
#define PASMIN_PFR 1.0E-3//Low resolution


#define MAX_NUM_TRACING 80000
#define NRELAX 3000

/* ========== Precision Settings ========== */
#define EPSILON_8 1.0E-8
#define EPSILON_10 1.0E-10
#define EPSILON_12 1.0E-12
#define EPSILON_15 1.0E-15

/* ========== 2DGRID and 3DGRID Settings ========== */
#define DEFAULT_POL_MARGIN 30
#define DEFAULT_RAD_MARGIN 30
#define DEFAULT_TOR_MARGIN 0


//for the inner and outer lege, the original resolution minus nlast/nfirst should larger than MIN_RESOLTUTION
#define MIN_RESOLTUTION 5

#define GRID_ALIGNMENT 32


#endif /* CONFIG_H */