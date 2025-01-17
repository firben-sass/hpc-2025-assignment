#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOLERANCE 1e-6

// Function to read binary file into a 3D array
double* read_binary_file(const char* filename, size_t items) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        perror("Failed to open file");
        return NULL;
    }

    double* data = (double*)malloc(items * sizeof(double));
    if (!data) {
        perror("Memory allocation failed");
        fclose(file);
        return NULL;
    }

    size_t read_items = fread(data, sizeof(double), items, file);
    if (read_items != items) {
        fprintf(stderr, "Error: Expected %lu items, but read %lu\n", items, read_items);
        free(data);
        fclose(file);
        return NULL;
    }

    fclose(file);
    return data;
}

// Function to compare two arrays
int compare_arrays(double* arr1, double* arr2, size_t items, double tolerance) {
    for (size_t i = 0; i < items; i++) {
        if (fabs(arr1[i] - arr2[i]) > tolerance) {
            printf("Difference at index %lu: %f vs %f\n", i, arr1[i], arr2[i]);
            return 0; // Not approximately equal
        }
    }
    return 1; // Arrays are approximately equal
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <file1> <file2>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char* file1 = argv[1];
    const char* file2 = argv[2];
    size_t num = 64; // Adjust based on your 3D array dimensions
    size_t items = num * num * num;

    double* data1 = read_binary_file(file1, items);
    double* data2 = read_binary_file(file2, items);

    if (!data1 || !data2) {
        free(data1);
        free(data2);
        return EXIT_FAILURE;
    }

    if (compare_arrays(data1, data2, items, TOLERANCE)) {
        printf("Test Passed: Outputs are approximately identical.\n");
    } else {
        printf("Test Failed: Outputs differ.\n");
    }

    free(data1);
    free(data2);
    return EXIT_SUCCESS;
}
