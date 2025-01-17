/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>

double *** jacobi(double ***u, double ***v, double ***f, int N, double *ki_to_xz, double *j_to_y) {
    // Using Jacobian method to solve the Poisson problem
    double factor = 1.0/6.0;
    double delta = 2.0/N;

    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            for (int k = 1; k < N+1; k++) {
                v[i][j][k] = factor * (u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta*delta * f[i][j][k]);
            }
        }
    }
    return v;
}
