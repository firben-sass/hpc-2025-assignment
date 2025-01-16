#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include "print.h"
#include "alloc3d.h"
#include "define_u_f.h"
#include "jacobi.h"
#include "gauss_seidel.h"


double*** deepcopy3DArray(double*** source, double *** copy, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                copy[i][j][k] = source[i][j][k];
            }
        }
    }
}

bool compare_3d_arrays(double *** arr1, double *** arr2, int N, double tol) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                if (fabs(arr1[i][j][k] - arr2[i][j][k]) > tol) {
                    return false;
                }
            }
        }
    }
    return true;
}

int test_compare_3d_arrays(double *** array_j, double *** array_gs, int N, double tol) {
    if (compare_3d_arrays(array_j, array_gs, N+2, tol)) {
        printf("Test Passed: Arrays are approximately identical.\n");
        return 1;
    } else {
        printf("Test Failed: Arrays are not approximately identical.\n");
        return 0;
    }
}

int main() {
    int 	N = 10;
    int 	iter_max = 10000;
    double	tolerance = 0.001;
    double	start_T = 0;
    double 	***u_0 = NULL;
    double  ***u_1 = NULL;
    double  ***out_gs  = NULL;
    double  ***out_j   = NULL;
    double  ***f   = NULL;


    // allocate memory
    if ( (u_0 = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u_0: allocation failed");
        exit(-1);
    }
    if ( (u_1 = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u_1: allocation failed");
        exit(-1);
    }
    if ( (f = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array f: allocation failed");
        exit(-1);
    }
    if ( (out_gs = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array out_gs: allocation failed");
        exit(-1);
    }
    if ( (out_j = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array out_j: allocation failed");
        exit(-1);
    }

    define_u(u_0, N);
    define_u(u_1, N);
    define_f(f, N);

    printf("Running Jacobi\n\n");
    seq_jacobi(u_0, u_1, f, N, iter_max);
    deepcopy3DArray(u_1, out_j, N);

    define_u(u_0, N);
    define_u(u_1, N);
    define_f(f, N);

    printf("Running Gauss Seidel\n\n");
    seq_gauss_seidel(u_0, u_1, f, N, iter_max);
    deepcopy3DArray(u_1, out_gs, N);

    // Print out_j array
    printf("out_j:\n");
    for (int i = 0; i < N+2; i++) {
        printf("2D Array at index %d:\n", i);
        for (int j = 0; j < N+2; j++) {
            for (int k = 0; k < N+2; k++) {
                printf("%#.6g ", out_j[i][j][k]); // Print each number with a minimum width of 3
            }
            printf("\n"); // Newline for rows
        }
        printf("\n"); // Newline between 2D arrays
    }

    // Print out_gs array
    printf("out_gs:\n");
    for (int i = 0; i < N+2; i++) {
        printf("2D Array at index %d:\n", i);
        for (int j = 0; j < N+2; j++) {
            for (int k = 0; k < N+2; k++) {
                printf("%#.6g ", out_gs[i][j][k]); // Print each number with a minimum width of 3
            }
            printf("\n"); // Newline for rows
        }
        printf("\n"); // Newline between 2D arrays
    }

    int result = test_compare_3d_arrays(out_j, out_gs, N, tolerance);
    
    // de-allocate memory
    free_3d(u_0);
    free_3d(u_1);
    free_3d(f);
    free_3d(out_gs);
    free_3d(out_j);

    return 0;
}
