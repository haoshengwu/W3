#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>

#define MAX_CONFIG_ITEMS 128
#define MAX_LINE_LENGTH 256
#define MAX_SECTION_NAME 128
#define MAX_KEY_NAME 256
#define MAX_VALUE_LENGTH 256
#define MAX_ZONES 16
#define MAX_XPOINT_NUMBER 4
//  Define a class which contain all the input parameters:
//    
//  equlibrium_file: the equlibrium file, which is assumed to be .geqdsk file.
//  
//  others will be added step by step.

typedef struct {
  char equilibrium_file[50];
  char topology_type[20];
  //The estimation of X-point position. Current it only consider ONE X point, will be updated.
  double xpt_estimation[2];  // unit is meter.

}  InputPara;

void init_inputpara(InputPara* input);
void print_inputpara(InputPara* input);


/*****************
*  INPUT MODULE  *
******************/

// ==================== Basic data structure ====================
typedef struct {
    char section[MAX_SECTION_NAME];
    char key[MAX_KEY_NAME];
    char value[MAX_VALUE_LENGTH];
} ConfigItem;

typedef struct {
    ConfigItem items[MAX_CONFIG_ITEMS];
    int count;
} ConfigData;


// ==================== Basic functions ====================

int parse_config_file(const char* filename, ConfigData* config);
char* trim_whitespace(char* str);
void print_config_data(const ConfigData* config);
void print_config_data_debug(const ConfigData* config);

// ==================== Parameter retrieval functions ====================
char* get_config_string(ConfigData* config, const char* section, 
                        const char* key, const char* default_val);
int get_config_int(ConfigData* config, const char* section, 
                   const char* key, int default_val);
double get_config_double(ConfigData* config, const char* section, 
                       const char* key, double default_val);
int get_config_bool(ConfigData* config, const char* section, 
                  const char* key, int default_val);

int has_config_key(ConfigData* config, const char* section, const char* key);

// ==================== Array parsing functions ====================
int parse_comma_separated_doubles(const char* str, double* array, int max_count);
int parse_comma_separated_ints(const char* str, int* array, int max_count);
int parse_zone_names(const char* str, char names[][32], int max_zones);



typedef struct {
    char geqdsk_file[256];
    char divgeo_trg_file[256];
} FileConfig;

typedef struct {
    char topology[32];
    int zone_number;
    char zone_names[MAX_ZONES][32];
    int xpoint_number; 
    double xpoint_r_est[MAX_XPOINT_NUMBER];
    double xpoint_z_est[MAX_XPOINT_NUMBER];
    double xpoint_delta;
    double guard_len_head[MAX_ZONES];
    double guard_len_tail[MAX_ZONES];
    double pasmin[MAX_ZONES];
    int neu_expand_number[MAX_ZONES];

} Grid2DConfig;

typedef struct {
    char toroidal_file[256];
    int toroidal_cell_number;
    double toroidal_delta;
} Grid3DConfig;

typedef struct {
    FileConfig file_config;
    Grid2DConfig grid2d_config;
    Grid3DConfig grid3d_config;
} W3Config;


// Load W3 configuration from file
int load_w3_config(const char* filename, W3Config* config);
// Validate W3 configuration
int validate_w3_config(const W3Config* config);
// Print W3 configuration
void print_w3_config(const W3Config* config);
// Free W3 configuration (currently no dynamic allocation, but good practice)
void free_w3_config(W3Config* config);

#endif

