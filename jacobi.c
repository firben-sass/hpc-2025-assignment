/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdio.h>
#include <omp.h>


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
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            for (int k = 1; k < N+1; k++) {
                v[i][j][k] = factor * (u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + delta*delta * f[i][j][k]);
            }
        }
    }
    return v;
}


// int jacobi(double *** u_0, double *** u_1, double *** f, double Delta, int N, int P)
// {
//     for (int p = 0; p < P; p++)
//     {
//         for (int i = 1; i < N+1; i++)
//         {
//             for (int j = 1; j < N+1; j++)
//             {
//                 for (int k = 1; k < N+1; k++)
//                 {
//                     u_1[i][j][k] = (u_0[i-1][j][k] + u_0[i+1][j][k] + u_0[i][j-1][k] + u_0[i][j+1][k] + u_0[i][j][k-1] + u_0[i][j][k+1] + Delta*Delta*f[i][j][k]) / 6;
//                 }
//             }
//         }

//         u_0 = u_1;
//     }

//     return 0;
// }
