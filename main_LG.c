/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#define N_DEFAULT 100

int
main(int argc, char *argv[]) {

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
        perror("array u: allocation failed");
        exit(-1);
    }

    if ( (f_matrix = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array f: allocation failed");
        exit(-1);
    }
    // define copy of u for Jacobi
    if ( (u_old = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array v: allocation failed");
        exit(-1);
    }

    // Create helper vector for the pointer gymnastics 
    double *ki_to_xz = (double *) malloc((N+1) * sizeof(double));
    double *j_to_y = (double *) malloc((N+1) * sizeof(double));
    double step_size = 2.0 / N;
    double x, y, z;
    ki_to_xz[1] = -1;
    j_to_y[1] = 1;

    for (int i = 2; i < N+1; i++) {
        ki_to_xz[i] = ki_to_xz[i-1] + step_size;
        j_to_y[i] = j_to_y[i-1] - step_size;
    }
    // Maybe add a small sanity check here that the end points are indeed -1 and 1

    // Loop through and initialize f and the boundary values of u
    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < N+2; j++) {
            for (int k = 0; k < N+2; k++) {
                if (i == 0 || i == N+1 || j == 0 || k == 0 || k == N+1) { // Can probably be optimized!
                    u[i][j][k] = 20.0;
                } else if (j == N+1) {
                    u[i][j][k] = 0.0;
                } else {
                    u[i][j][k] = start_T;
                }
                f_matrix[i][j][k] = f(ki_to_xz[k], j_to_y[j], ki_to_xz[i]);
            }
        }
    }

    // Iterations
    int iter = 0;
    double diff = 0.0;
    while (diff > tolerance && iter < iter_max) {
        // Copy u to u_old
        u_old = memcpy(u_old, u, (N+2)*(N+2)*(N+2)*sizeof(double));

        // JACOBI
        u = jacobi(u, v, f_matrix, N, ki_to_xz, j_to_y);
        // Calculate the difference between the new and the old u
        diff = 0.0;
        for (int i = 1; i < N+1; i++) {
            for (int j = 1; j < N+1; j++) {
                for (int k = 1; k < N+1; k++) {
                    // Compute difference using 2-norm
                    diff += (u[i][j][k] - u_old[i][j][k]) * (u[i][j][k] - u_old[i][j][k]);
                }
            }
        }
        iter++;
    }

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
    free_3d(u);
    free_3d(f_matrix);

    return(0);
}

int f(double x, double y, double z) {
    if (x >= -1 && x <= -0.375 && y >= -1 && y <= -0.5 && z >=-2/3 && z <= 0) {
        return 200;
    } else {
        return 0;
    }
}
