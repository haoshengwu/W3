#include "utils.h"

#include <stdlib.h>
#include <stdio.h>

void write_array(const double* arr, int n, char* filename)
{
    if (!arr || n <= 0 || !filename)
    {
        fprintf(stderr, "Invalid input to write_array.\n");
        exit(EXIT_FAILURE);
    }

    FILE* fp = fopen(filename, "w");
    if (!fp)
    {
        perror("Failed to open file for writing");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; ++i)
    {
        fprintf(fp, "%.15e\n", arr[i]); 
    }

    fclose(fp);
    printf("Finish write array to %s\n", filename);
}

void write_array2f(const double* x, const double* y, int n, char* filename)
{
    if (!x || !y || n <= 0 || !filename)
    {
        fprintf(stderr, "Invalid input to write_array2f.\n");
        exit(EXIT_FAILURE);
    }

    FILE* fp = fopen(filename, "w");
    if (!fp)
    {
        perror("Failed to open file in write_array2f");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; ++i)
    {
        fprintf(fp, "%.15e %.15e\n", x[i], y[i]);
    }

    fclose(fp);
    printf("Finish write array2f to %s\n", filename);
}