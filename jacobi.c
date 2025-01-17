/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdio.h>

double *** jacobi(double ***u, double ***v, double ***f, int N) {
    // Using Jacobian method to solve the Poisson problem
    double factor = 1.0/6.0;
    double delta = 2.0/( (double) (N + 1));

    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            for (int k = 1; k < N+1; k++) {
                v[i][j][k] = factor * (u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta*delta * f[i][j][k]);
            }
        }
    }
    return v;
}

double *** jacobi_par(double ***u, double ***v, double ***f, int N) {
    // Using Jacobian method to solve the Poisson problem
    double factor = 1.0/6.0;
    double delta = 2.0/( (double) (N + 1));
    // parallelization
    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            for (int k = 1; k < N+1; k++) {
                v[i][j][k] = factor * (u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta*delta * f[i][j][k]);
            }
        }
    }
    return v;
}