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

// Define source function f (radiator)
int f(double x, double y, double z) {
    if (x >= -1 && x <= -0.375 && y >= -1 && y <= -0.5 && z >=-2/3 && z <= 0) {
        return 200;
    } else {
        return 0;
    }
}

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
    double start_time, end_time;

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

    // Create helper vector for the pointer gymnastics 
    double *x = (double *) malloc((N+1) * sizeof(double));
    double *y = (double *) malloc((N+1) * sizeof(double));
    double *z = (double *) malloc((N+1) * sizeof(double));
    double step_size = 2.0 / (N-1);
    x[1] = -1;
    y[1] = 1;
    z[1] = -1;

    for (int i = 2; i < N+1; i++) {
        x[i] = x[i-1] + step_size;
        y[i] = y[i-1] - step_size;
        z[i] = z[i-1] + step_size;
    }

    // Loop through and initialize f and the boundary values of u (and initializing v = u)
    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < N+2; j++) {
            for (int k = 0; k < N+2; k++) {
                if (i == 0 || i == N+1 || j == 0 || k == 0 || k == N+1) { // Can probably be optimized!
                    u[i][j][k] = 20.0;
                    v[i][j][k] = 20.0;
                } else if (j == N+1) {
                    u[i][j][k] = 0.0;
                    v[i][j][k] = 0.0;
                } else {
                    u[i][j][k] = start_T;
                    v[i][j][k] = start_T;
                }
                f_matrix[i][j][k] = f(x[k], y[j], z[i]);
            }
        }
    }

    // Iterations
    int iter = 0;
    double diff = 9999999.0;
    double diff_1;
    // Running the iterations until the difference is less than the tolerance or max number of iterations is reached
    start_time = omp_get_wtime();
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
        v = jacobi(u, v, f_matrix, N);

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
                    // printf("diff[%d][%d][%d] = %f\n", i, j, k, (u[i][j][k] - u_old[i][j][k]));
                }
            }
        }
        diff = diff_1;
        // printf("Iteration: %d, diff: %f\n", iter, diff);
        iter++;
    }
    end_time = omp_get_wtime();
    printf("Used %d iterations, diff = %.6f, time taken: %.8f\n", iter, diff, end_time - start_time);
    // printf("Time taken: %f\n", end_time - start_time);

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
    // printf("Deallocating memory\n");
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

    free(x); free(y); free(z);

    return(0);
}


