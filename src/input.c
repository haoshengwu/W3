
#include "input.h"

// Here are the default values of the input parameters for building the mesh
// At the begining, we assume is single null. It the future, it should suport
// different configurations.


// Current use DTT values, will be updated in the future.
void init_inputpara(InputPara* input){

  strcpy(input -> equilibrium_file,"example.geqdsk");
  strcpy(input -> topology_type,"sinlge_null");
  //for DTT
  input->xpt_estimation[0] = 2.16;  //  R coordinate
  input->xpt_estimation[1] = -1.11; //  Z coordinate
  //for SPARC
  //input->xpt_estimation[0] = 1.5417;  //  R coordinate
  //input->xpt_estimation[1] = -1.1167; //  Z coordinate
}
void print_inputpara(InputPara* input){

  printf("Equilibirum file name: %s\n", input -> equilibrium_file);
  printf("The topology of equilibrium file is: %s\n", input -> topology_type);
  
}

//====================================================================================

static char* string_duplicate(const char* str) 
{
    if (!str) return NULL;
    
    size_t len = strlen(str) + 1;  // +1 for null terminator
    char* copy = malloc(len);
    
    if (!copy) {
        return NULL;  // Memory allocation failed
    }

    strcpy(copy, str);
    return copy;
}


int parse_config_file(const char* filename, ConfigData* config) {
   // 1. Open file
   FILE* fp = fopen(filename, "r");
   if (!fp) {
       fprintf(stderr, "Cannot open config file: %s\n", filename);
       return -1;  // File open failed
   }
   
   // 2. Initialize variables
   char line[MAX_LINE_LENGTH];           // Buffer to store each line content
   char current_section[MAX_SECTION_NAME] = "";  // Current section name
   config->count = 0;                    // Reset parsed config item count
   
   // 3. Read file line by line
   while (fgets(line, sizeof(line), fp)) {
       
       // 4. Handle comments - remove content after #
       char* comment = strchr(line, '#');     // Find first # character
       if (comment) *comment = '\0';          // Truncate string
       
       // 5. Remove leading and trailing whitespace
       char* trimmed = trim_whitespace(line);
       
       // 6. Skip empty lines
       if (strlen(trimmed) == 0) continue;
       
       // 7. Handle section line [SECTION_NAME]
       if (trimmed[0] == '[') {
           char* end = strchr(trimmed, ']');  // Find closing ]
           if (end) {
               *end = '\0';                   // Remove ]
               strcpy(current_section, trimmed + 1);  // Copy section name (skip opening [)
               trim_whitespace(current_section);      // Remove whitespace from section name
           }
           continue;  // Done with section, continue to next line
       }
       
       // 8. Handle key=value line
       char* equals = strchr(trimmed, '=');   // Find = sign
       if (equals && config->count < MAX_CONFIG_ITEMS) {
           *equals = '\0';                    // Split at =, separate key and value
           char* key = trim_whitespace(trimmed);       // Get key and trim whitespace
           char* value = trim_whitespace(equals + 1);  // Get value and trim whitespace
           
           // 9. Store into ConfigData structure
           ConfigItem* item = &config->items[config->count];
           strcpy(item->section, current_section);    // Copy current section name
           strcpy(item->key, key);                   // Copy key
           strcpy(item->value, value);               // Copy value
           config->count++;                          // Increment count
       }
   }
   
   // 10. Close file and return success
   fclose(fp);
   return 0;
}

char* trim_whitespace(char* str) {
   // Check for null pointer
   if (!str) return NULL;
   
   // Remove leading whitespace characters
   // Move pointer forward while current character is whitespace
   while (isspace((unsigned char)*str)) str++;
   
   // If string contains only whitespace characters
   if (*str == 0) return str;
   
   // Remove trailing whitespace characters
   // Find the last character of the string
   char* end = str + strlen(str) - 1;
   
   // Move backwards from end while character is whitespace
   // Stop when we reach the beginning or find non-whitespace
   while (end > str && isspace((unsigned char)*end)) end--;
   
   // Null-terminate the string after the last non-whitespace character
   *(end + 1) = '\0';
   
   // Return pointer to the trimmed string
   return str;
}

void print_config_data(const ConfigData* config) {
   if (!config) {
       printf("ConfigData is NULL\n");
       return;
   }
   
   printf("=== Configuration Data ===\n");
   printf("Total items: %d\n\n", config->count);
   
   if (config->count == 0) {
       printf("No configuration items found.\n");
       return;
   }
   
   // Print items grouped by section
   char current_section[MAX_SECTION_NAME] = "";
   
   for (int i = 0; i < config->count; i++) {
       const ConfigItem* item = &config->items[i];
       
       // Print section header when section changes
       if (strcmp(current_section, item->section) != 0) {
           strcpy(current_section, item->section);
           
           // Add blank line before new section (except for first section)
           if (i > 0) printf("\n");
           
           if (strlen(item->section) > 0) {
               printf("[%s]\n", item->section);
           } else {
               printf("[Global]\n");  // For items without section
           }
       }
       
       // Print key-value pair with proper formatting
       printf("  %-20s = %s\n", item->key, item->value);
   }
   
   printf("\n=== End of Configuration ===\n");
}

void print_config_data_debug(const ConfigData* config) {
    if (!config) {
        printf("ConfigData is NULL\n");
        return;
    }

    printf("=== Configuration Data (Debug Mode) ===\n");
    printf("Structure address: %p\n", (void*)config);
    printf("Total items: %d\n", config->count);
    printf("Max capacity: %d\n\n", MAX_CONFIG_ITEMS);
    
    if (config->count == 0) {
        printf("No configuration items found.\n");
        return;
    }
    
    // Print each item with index and detailed info
    for (int i = 0; i < config->count; i++) {
        const ConfigItem* item = &config->items[i];
        
        printf("Item[%2d]:\n", i);
        printf("  Section : '%s' (len=%zu)\n", item->section, strlen(item->section));
        printf("  Key     : '%s' (len=%zu)\n", item->key, strlen(item->key));
        printf("  Value   : '%s' (len=%zu)\n", item->value, strlen(item->value));
        printf("  Address : %p\n", (void*)item);
        printf("  ---\n");
    }
    
    printf("=== End of Debug Info ===\n");
}

// Get string parameter from configuration
char* get_config_string(ConfigData* config, const char* section, 
                      const char* key, const char* default_val) {
   // Check input parameters
   if (!config || !section || !key) {
       return (char*)default_val;
   }
   
   // Search through all config items for matching section and key
   for (int i = 0; i < config->count; i++) {
       if (strcmp(config->items[i].section, section) == 0 &&
           strcmp(config->items[i].key, key) == 0) {
           return config->items[i].value;  // Return found value
       }
   }
   
   // Return default if not found
   return (char*)default_val;
}

// Get integer parameter with validation
int get_config_int(ConfigData* config, const char* section, 
                 const char* key, int default_val) {
   char* str = get_config_string(config, section, key, NULL);
   
   if (!str) {
       return default_val;
   }
   
   // Use strtol for safe conversion
   char* endptr;
   errno = 0;  // Reset errno
   long val = strtol(str, &endptr, 10);
   
   // Check if conversion was successful
   if (endptr == str || *endptr != '\0') {
       fprintf(stderr, "Warning: '%s' is not a valid integer for [%s].%s, using default %d\n", 
               str, section, key, default_val);
       return default_val;
   }
   
   // Check for overflow
   if (errno == ERANGE || val > INT_MAX || val < INT_MIN) {
       fprintf(stderr, "Warning: Value '%s' out of range for [%s].%s, using default %d\n", 
               str, section, key, default_val);
       return default_val;
   }
   
   return (int)val;
}

// Get double parameter with validation
double get_config_double(ConfigData* config, const char* section, 
                       const char* key, double default_val) {
   char* str = get_config_string(config, section, key, NULL);
   
   if (!str) {
       return default_val;
   }
   
   // Use strtod for safe conversion
   char* endptr;
   errno = 0;  // Reset errno
   double val = strtod(str, &endptr);
   
   // Check if conversion was successful
   if (endptr == str || *endptr != '\0') {
       fprintf(stderr, "Warning: '%s' is not a valid number for [%s].%s, using default %.6g\n", 
               str, section, key, default_val);
       return default_val;
   }
   
   // Check for overflow or underflow
   if (errno == ERANGE) {
       fprintf(stderr, "Warning: Value '%s' out of range for [%s].%s, using default %.6g\n", 
               str, section, key, default_val);
       return default_val;
   }
   
   return val;
}

// Get boolean parameter with validation
int get_config_bool(ConfigData* config, const char* section, 
                  const char* key, int default_val) {
   char* str = get_config_string(config, section, key, NULL);
   
   if (!str) {
       return default_val;
   }
   
   // Create lowercase copy for comparison
   char lower[64];
   int i = 0;
   
   // Convert to lowercase, handle longer strings
   while (str[i] && i < 63) {
       lower[i] = tolower((unsigned char)str[i]);
       i++;
   }
   lower[i] = '\0';
   
   // Check for true values
   if (strcmp(lower, "true") == 0 || 
       strcmp(lower, "yes") == 0 || 
       strcmp(lower, "1") == 0 || 
       strcmp(lower, "on") == 0 ||
       strcmp(lower, "enable") == 0 ||
       strcmp(lower, "enabled") == 0) {
       return 1;
   }
   
   // Check for false values
   if (strcmp(lower, "false") == 0 || 
       strcmp(lower, "no") == 0 || 
       strcmp(lower, "0") == 0 || 
       strcmp(lower, "off") == 0 ||
       strcmp(lower, "disable") == 0 ||
       strcmp(lower, "disabled") == 0) {
       return 0;
   }
   
   // Unrecognized value, give warning
   fprintf(stderr, "Warning: '%s' is not a valid boolean for [%s].%s, using default %s\n", 
           str, section, key, default_val ? "true" : "false");
   return default_val;
}

// Check if a key exists in configuration
int has_config_key(ConfigData* config, const char* section, const char* key) {
   // Check input parameters
   if (!config || !section || !key) {
       return 0;
   }
   
   // Search through all config items for matching entry
   for (int i = 0; i < config->count; i++) {
       if (strcmp(config->items[i].section, section) == 0 &&
           strcmp(config->items[i].key, key) == 0) {
           return 1;  // Found
       }
   }
   
   return 0;  // Not found
}

// ==================== Array parsing functions ====================

// Parse comma-separated double values
int parse_comma_separated_doubles(const char* str, double* array, int max_count) {
   // Check input parameters
   if (!str || !array || max_count <= 0) {
       return 0;
   }
   
   // Skip if string is empty
   if (strlen(str) == 0) {
       return 0;
   }
   
   // Create a copy of the string for tokenization
   char* copy = string_duplicate(str);
   if (!copy) {
       fprintf(stderr, "Error: Memory allocation failed in parse_comma_separated_doubles\n");
       return 0;
   }
   
   char* token = strtok(copy, ",");
   int count = 0;
   
   // Process each token
   while (token && count < max_count) {
       // Trim whitespace from token
       char* trimmed = trim_whitespace(token);
       
       // Convert to double with validation
       char* endptr;
       errno = 0;
       double val = strtod(trimmed, &endptr);
       
       // Check if conversion was successful
       if (endptr == trimmed || *endptr != '\0') {
           fprintf(stderr, "Warning: '%s' is not a valid number, skipping\n", trimmed);
       } else if (errno == ERANGE) {
           fprintf(stderr, "Warning: '%s' is out of range, skipping\n", trimmed);
       } else {
           array[count] = val;
           count++;
       }
       
       token = strtok(NULL, ",");
   }
   
   free(copy);
   return count;
}

// Parse comma-separated integer values
int parse_comma_separated_ints(const char* str, int* array, int max_count) {
   // Check input parameters
   if (!str || !array || max_count <= 0) {
       return 0;
   }
   
   // Skip if string is empty
   if (strlen(str) == 0) {
       return 0;
   }
   
   // Create a copy of the string for tokenization
   char* copy = string_duplicate(str);
   if (!copy) {
       fprintf(stderr, "Error: Memory allocation failed in parse_comma_separated_ints\n");
       return 0;
   }
   
   char* token = strtok(copy, ",");
   int count = 0;
   
   // Process each token
   while (token && count < max_count) {
       // Trim whitespace from token
       char* trimmed = trim_whitespace(token);
       
       // Convert to integer with validation
       char* endptr;
       errno = 0;
       long val = strtol(trimmed, &endptr, 10);
       
       // Check if conversion was successful
       if (endptr == trimmed || *endptr != '\0') {
           fprintf(stderr, "Warning: '%s' is not a valid integer, skipping\n", trimmed);
       } else if (errno == ERANGE || val > INT_MAX || val < INT_MIN) {
           fprintf(stderr, "Warning: '%s' is out of range, skipping\n", trimmed);
       } else {
           array[count] = (int)val;
           count++;
       }
       
       token = strtok(NULL, ",");
   }
   
   free(copy);
   return count;
}

// Parse comma-separated zone names
int parse_zone_names(const char* str, char names[][32], int max_zones) {
   // Check input parameters
   if (!str || !names || max_zones <= 0) {
       return 0;
   }
   
   // Skip if string is empty
   if (strlen(str) == 0) {
       return 0;
   }
   
   // Create a copy of the string for tokenization
   char* copy = string_duplicate(str);
   if (!copy) {
       fprintf(stderr, "Error: Memory allocation failed in parse_zone_names\n");
       return 0;
   }
   
   char* token = strtok(copy, ",");
   int count = 0;
   
   // Process each token
   while (token && count < max_zones) {
       // Trim whitespace from token
       char* trimmed = trim_whitespace(token);
       
       // Check if name is not empty after trimming
       if (strlen(trimmed) > 0) {
           // Copy name ensuring null termination
           strncpy(names[count], trimmed, 31);
           names[count][31] = '\0';  // Ensure null termination
           count++;
       } else {
           fprintf(stderr, "Warning: Empty zone name found, skipping\n");
       }

       token = strtok(NULL, ",");
   }
   
   free(copy);
   return count;
}

int load_w3_config(const char* filename, W3Config* config) {
    ConfigData raw_config = {0};
    
    // Parse configuration file
    if (parse_config_file(filename, &raw_config) != 0) 
    {
        fprintf(stderr, "Failed to parse config file: %s\n", filename);
        return -1;
    }
    
    // Initialize all arrays to zero
    memset(config, 0, sizeof(W3Config));
    
    // ==================== Load FILE section ====================
    strcpy(config->file_config.geqdsk_file, 
           get_config_string(&raw_config, "FILE", "GEDQSK_FILE", "example.gedqsk"));
    strcpy(config->file_config.divgeo_trg_file,
           get_config_string(&raw_config, "FILE", "DIVGEO_TRG_FILE", "example.trg"));
    
    // ==================== Load 2DGRID_CONTROL section ====================
    strcpy(config->grid2d_config.topology,
           get_config_string(&raw_config, "2DGRID_CONTROL", "TOPOLOGY", "SNL"));
    
    config->grid2d_config.zone_number = 
        get_config_int(&raw_config, "2DGRID_CONTROL", "ZONE_NUMBER", 3);
    
    // Parse zone names
    char* zone_name_str = get_config_string(&raw_config, "2DGRID_CONTROL", "ZONE_NAME", "SOL,PFR,CORE");
    parse_zone_names(zone_name_str, config->grid2d_config.zone_names, MAX_ZONES);
    
    // Load X-point configuration
    config->grid2d_config.xpoint_number = 
    get_config_int(&raw_config, "2DGRID_CONTROL", "XPOINT_NUMBER", 1);

    // Handle X-point R positions with strict validation
    char* xpoint_r_str = get_config_string(&raw_config, "2DGRID_CONTROL", "XPOINT_R_EST", "");
    if (strlen(xpoint_r_str) > 0) 
    {
        int r_count = parse_comma_separated_doubles(xpoint_r_str, 
                                                   config->grid2d_config.xpoint_r_est, 
                                                   MAX_XPOINT_NUMBER);
            // Check if number of R values matches XPOINT_NUMBER
        if (r_count != config->grid2d_config.xpoint_number) 
        {
            fprintf(stderr, "Error: XPOINT_R_EST has %d values but XPOINT_NUMBER is %d. "
                    "Number of values must match exactly.\n", 
                    r_count, config->grid2d_config.xpoint_number);
            exit(EXIT_FAILURE);  // Exit with error
        }
    } 
    else 
    {
        // If XPOINT_R_EST is not provided but XPOINT_NUMBER > 0
        if (config->grid2d_config.xpoint_number > 0) {
            fprintf(stderr, "Error: XPOINT_R_EST is required when XPOINT_NUMBER = %d\n", 
                config->grid2d_config.xpoint_number);
            exit(EXIT_FAILURE);  // Exit with error
        }
    }

// Handle X-point Z positions with strict validation
    char* xpoint_z_str = get_config_string(&raw_config, "2DGRID_CONTROL", "XPOINT_Z_EST", "");
    if (strlen(xpoint_z_str) > 0) 
    {
        int z_count = parse_comma_separated_doubles(xpoint_z_str, 
                                                    config->grid2d_config.xpoint_z_est, 
                                                    MAX_XPOINT_NUMBER);
        // Check if number of Z values matches XPOINT_NUMBER
        if (z_count != config->grid2d_config.xpoint_number) 
        {
            fprintf(stderr, "Error: XPOINT_Z_EST has %d values but XPOINT_NUMBER is %d. "
                    "Number of values must match exactly.\n", 
                    z_count, config->grid2d_config.xpoint_number);
        exit(EXIT_FAILURE);  // Exit with error
        }
    } 
    else 
    {
      // If XPOINT_Z_EST is not provided but XPOINT_NUMBER > 0
        if (config->grid2d_config.xpoint_number > 0) 
        {
            fprintf(stderr, "Error: XPOINT_Z_EST is required when XPOINT_NUMBER = %d\n", 
                   config->grid2d_config.xpoint_number);
            exit(EXIT_FAILURE);// Exit with error
        }
    }
    
    config->grid2d_config.xpoint_delta = 
        get_config_double(&raw_config, "2DGRID_CONTROL", "XPOINT_DELTA", 0.1);
    
    // Parse array parameters
    char* guard_head_str = get_config_string(&raw_config, "2DGRID_CONTROL", "GUARD_LEN_HEAD", "0.1,0.1,0.0");
    parse_comma_separated_doubles(guard_head_str, config->grid2d_config.guard_len_head, MAX_ZONES);
    
    char* guard_tail_str = get_config_string(&raw_config, "2DGRID_CONTROL", "GUARD_LEN_TAIL", "0.1,0.1,0.0");
    parse_comma_separated_doubles(guard_tail_str, config->grid2d_config.guard_len_tail, MAX_ZONES);
    
    char* pasmin_str = get_config_string(&raw_config, "2DGRID_CONTROL", "PASMIN", "1.0E-3,1.0E-3,1.0E-3");
    parse_comma_separated_doubles(pasmin_str, config->grid2d_config.pasmin, MAX_ZONES);

    char* neu_expand_number = get_config_string(&raw_config, "2DGRID_CONTROL", "NEU_EXPAND_NUMBER", "2,2");
    parse_comma_separated_ints(neu_expand_number, config->grid2d_config.neu_expand_number, MAX_ZONES);
    
    // ==================== Load 3DGRID_CONTROL section ====================
    strcpy(config->grid3d_config.toroidal_file,
           get_config_string(&raw_config, "3DGRID_CONTROL", "TOROIDAL_FILE", ""));
    
    config->grid3d_config.toroidal_cell_number = 
        get_config_int(&raw_config, "3DGRID_CONTROL", "TOROIDAL_CELL_NUMBER", 10);

    config->grid3d_config.toroidal_delta = 
        get_config_double(&raw_config, "3DGRID_CONTROL", "TOROIDAL_DELTA", 1.0);
    
    return 0;
}


// Validate physical parameters for each zone
int validate_w3_config(const W3Config* config)
{
  //TODO
  return 0;
}

void print_w3_config(const W3Config* config) {
    if (!config) {
        printf("W3Config is NULL\n");
        return;
    }


    printf("======================================\n");
    printf("=          W3 CONFIGURATION          =\n");
    printf("======================================\n\n");
    printf("\n");
    
    // Print FILE section
    printf("[FILE]\n");
    printf("  GEDQSK_FILE     : %s\n", config->file_config.geqdsk_file);
    printf("  DIVGEO_TRG_FILE : %s\n", config->file_config.divgeo_trg_file);
    printf("\n");
    
    // Print 2DGRID_CONTROL section
    printf("[2DGRID_CONTROL]\n");
    printf("  TOPOLOGY           : %s\n", config->grid2d_config.topology);
    printf("  ZONE_NUMBER        : %d\n", config->grid2d_config.zone_number);
    
    printf("  ZONE_NAMES         : ");
    for (int i = 0; i < config->grid2d_config.zone_number; i++) {
        printf("%s", config->grid2d_config.zone_names[i]);
        if (i < config->grid2d_config.zone_number - 1) printf(", ");
    }
    printf("\n");
    
    // Print X-point information
    printf("  XPOINT_NUMBER      : %d\n", config->grid2d_config.xpoint_number);
    printf("  XPOINT_R_EST       : ");
    for (int i = 0; i < config->grid2d_config.xpoint_number; i++) {
        printf("%.3f", config->grid2d_config.xpoint_r_est[i]);
        if (i < config->grid2d_config.xpoint_number - 1) printf(", ");
    }
    printf("\n");
    
    printf("  XPOINT_Z_EST       : ");
    for (int i = 0; i < config->grid2d_config.xpoint_number; i++) {
        printf("%.3f", config->grid2d_config.xpoint_z_est[i]);
        if (i < config->grid2d_config.xpoint_number - 1) printf(", ");
    }
    printf("\n");
    
    printf("  XPOINT_DELTA       : %.6f\n", config->grid2d_config.xpoint_delta);
    
    // Print array parameters
    printf("  GUARD_LEN_HEAD     : ");
    for (int i = 0; i < config->grid2d_config.zone_number; i++) {
        printf("%.3f", config->grid2d_config.guard_len_head[i]);
        if (i < config->grid2d_config.zone_number - 1) printf(", ");
    }
    printf("\n");
    
    printf("  GUARD_LEN_TAIL     : ");
    for (int i = 0; i < config->grid2d_config.zone_number; i++) {
        printf("%.3f", config->grid2d_config.guard_len_tail[i]);
        if (i < config->grid2d_config.zone_number - 1) printf(", ");
    }
    printf("\n");
    
    printf("  PASMIN             : ");
    for (int i = 0; i < config->grid2d_config.zone_number; i++) {
        printf("%.3e", config->grid2d_config.pasmin[i]);
        if (i < config->grid2d_config.zone_number - 1) printf(", ");
    }
    printf("\n");

    printf("  NEU_EXPAND_NUMBER  : ");
    
    for (int i = 0; i < config->grid2d_config.zone_number; i++) {
        printf("%d", config->grid2d_config.neu_expand_number[i]);
        if (i < config->grid2d_config.zone_number - 1) printf(", ");
    }
    
    printf("\n\n");
    
    // Print 3DGRID_CONTROL section
    printf("[3DGRID_CONTROL]\n");
    printf("  TOROIDAL_FILE          : %s\n", config->grid3d_config.toroidal_file);
    printf("  TOROIDAL_CELL_NUMBER   : %d\n", config->grid3d_config.toroidal_cell_number);
    printf("  TOROIDAL_DELTA         : %.2f\n", config->grid3d_config.toroidal_delta);

    printf("\n");
}

// Free W3 configuration (currently no dynamic allocation, but good practice)
void free_w3_config(W3Config* config) 
{
    if (!config) return;    
    // Currently no dynamic memory to free, but this function
    // provides a place for future cleanup if needed
    memset(config, 0, sizeof(W3Config));
}