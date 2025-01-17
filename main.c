/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc3d.h"
#include "print.h"
#include <omp.h>

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#define N_DEFAULT 100

void print_3d_matrix(double ***matrix, int N) {
    for (int i = 0; i < N+2; i++) {
        printf("i = %d\n", i);
        for (int j = 0; j < N+2; j++) {
            printf("[%d][%d] = ", i, j);
            for (int k = 0; k < N+2; k++) {
                printf(" %f ", matrix[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");    
}

void print_matrix(double ***matrix, int N) {
// Print u array
    for (int k = 0; k < N+2; k++) {
        printf("2D Array at x[%d]:\n", k);
        for (int j = 0; j < N+2; j++) {
            for (int i = 0; i < N+2; i++) {
                printf("%#.6g ", matrix[i][j][k]); // Print each number with a minimum width of 3
            }
            printf("\n"); // Newline for rows
        }
        printf("\n"); // Newline between 2D arrays
    }
}

int main(int argc, char *argv[]) {

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;
    double  ***f_matrix = NULL;
    double  ***u_old = NULL;
    double  ***v = NULL;
    double  start_time, end_time;
    double  x,y,z;

    /* get the parameters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

    // allocate memory
    if ( (u = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    if ( (v = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array v: allocation failed");
        exit(-1);
    }

    if ( (f_matrix = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array f_matrix: allocation failed");
        exit(-1);
    }

    // define copy of u to u_old for computing difference
    if ( (u_old = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u_old: allocation failed");
        exit(-1);
    }
    
    double step_size = 2.0 / (N+1);

    // Loop through and initialize f and the boundary values of u (and initializing v = u)
    for (int i = 0; i < N+2; i++) {
        z = -1 + i*step_size;
        for (int j = 0; j < N+2; j++) {
            y = 1 - j*step_size;
            for (int k = 0; k < N+2; k++) {
                x = -1 + k*step_size;
                // Fill in f-matrix
                if (x >= -1.0 && x <= -3.0 / 8.0 && y >= -1.0 && y <= -0.5 &&
                    z >= -2.0 / 3.0 && z <= 0.0) {
                    f_matrix[i][j][k] = 200.0; // Radiator region
                } else {
                    f_matrix[i][j][k] = 0.0;   // Outside radiator region
                }
                // Fill in u and v
                if (x == -1.0 || x == 1.0 || z == -1.0 ||z == 1.0 || y == 1.0) {
                    u[i][j][k] = 20.0;
                    v[i][j][k] = 20.0;
                } else if (y == -1.0) {
                    u[i][j][k] = 0.0;
                    v[i][j][k] = 0.0;
                } else {
                    u[i][j][k] = (double) start_T;
                    v[i][j][k] = (double) start_T;
                }
            }
        }
    }
    // printf("f after initialization: ");
    // print_matrix(f_matrix, N);

    // Iterations
    int iter = 0;
    double diff = 9999999.0;
    double diff_1;
    // Running the iterations until the difference is less than the tolerance or max number of iterations is reached
    double total_time = 0.0;
    while (diff > tolerance && iter < iter_max) {
        // Copy u to u_old
        // Allocate a temporary 3D array
        double ***temp = malloc_3d(N+2, N+2, N+2);
        if (temp == NULL) {
            perror("Temp array allocation failed");
            exit(-1);
        }

        // Copy the content of u to temp
        for (int i = 0; i < N+2; i++) {
            for (int j = 0; j < N+2; j++) {
                for (int k = 0; k < N+2; k++) {
                    temp[i][j][k] = u[i][j][k];
                }
            }
        }

        // Copy the content of temp to u_old
        for (int i = 0; i < N+2; i++) {
            for (int j = 0; j < N+2; j++) {
                for (int k = 0; k < N+2; k++) {
                    u_old[i][j][k] = temp[i][j][k];
                }
            }
        }

        // Free the temporary array
        free_3d(temp);

        // Perform the Jacobi or Gauss-Seidel iteration
        #ifdef _JACOBI
        {
        start_time = omp_get_wtime();
        v = jacobi_par(u, v, f_matrix, N);
        end_time = omp_get_wtime();
        total_time += end_time - start_time;
        }

        // Copy the content of v back to u
        double ***tmp = malloc_3d(N+2, N+2, N+2);
        if (tmp == NULL) {
            perror("Tmp array allocation failed");
            exit(-1);
        }

        // Copy the content of v to tmp
        for (int i = 0; i < N+2; i++) {
            for (int j = 0; j < N+2; j++) {
                for (int k = 0; k < N+2; k++) {
                    tmp[i][j][k] = v[i][j][k];
                }
            }
        }

        // Copy the content of tmp to u
        for (int i = 0; i < N+2; i++) {
            for (int j = 0; j < N+2; j++) {
                for (int k = 0; k < N+2; k++) {
                    u[i][j][k] = tmp[i][j][k];
                }
            }
        }

        // Free the temporary array
        free_3d(tmp);
        #endif

        #ifdef _GAUSS_SEIDEL
        u = gauss_seidel(u, f_matrix, N);
        #endif

        // Calculate the difference between the new and the old u
        diff_1 = 0.0;
        for (int i = 1; i < N+1; i++) {
            for (int j = 1; j < N+1; j++) {
                for (int k = 1; k < N+1; k++) {
                    // Compute difference using 2-norm
                    diff_1 += (u[i][j][k] - u_old[i][j][k]) * (u[i][j][k] - u_old[i][j][k]);
                }
            }
        }
        diff = diff_1;
        iter++;
    }
    printf("Used %d iterations, diff = %.6f, time taken: %.8f\n", iter, diff, total_time);
    // printf("Time taken: %f\n", end_time - start_time);

    // print_matrix(u, N);
    // print_3d_matrix(u, N);

    // dump  results if wanted 
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    if (u != NULL) {
        free_3d(u);
    }
    if (v != NULL) {
        free_3d(v);
    }
    if (u_old != NULL) {
        free_3d(u_old);
    }
    if (f_matrix != NULL) {
        free_3d(f_matrix);
    }

    return(0);
}