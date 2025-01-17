/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>
#include <omp.h>


void _seq_rec_gauss_seidel(double *** u_0, double *** u_1, double *** f, int N, int i, int j, int k)
{   
    // If in a wall
    if (i == 0 || i == N+1 || j == 0 || j == N+1 || k == 0 || k == N+1)
        return;

    _seq_rec_gauss_seidel(u_0, u_1, f, N, i-1, j, k);
    double u_1_1_0_0 = u_1[i-1][j][k];
    _seq_rec_gauss_seidel(u_0, u_1, f, N, i, j-1, k);
    double u_1_0_1_0 = u_1[i][j-1][k];
    _seq_rec_gauss_seidel(u_0, u_1, f, N, i, j, k-1);
    double u_1_0_0_1 = u_1[i][j][k-1];

    double u_0_1_0_0 = u_0[i+1][j][k];
    double u_0_0_1_0 = u_0[i][j+1][k];
    double u_0_0_0_1 = u_0[i][j][k+1];

    double delta = 2.0/N;
    u_1[i][j][k] = (u_1_1_0_0 + u_1_0_1_0 + u_1_0_0_1 + u_0_1_0_0 + u_0_0_1_0 + u_0_0_0_1 + delta*delta*f[i][j][k]) / 6.0;

    return;
}

void seq_rec_gauss_seidel(double *** u_0, double *** u_1, double *** f, int N, int P)
{
    for (int p; p < P; p++)
    {
        _seq_rec_gauss_seidel(u_0, u_1, f, N, N, N, N);
        u_0 = u_1;
    }

    return;
}

void seq_gauss_seidel(double ***u_1, double ***f, int N, int P)
{
    // Grid spacing
    double delta = 2.0 / N;

    // Perform P iterations
    for (int p = 0; p < P; p++) {
        // Iterate over all grid points, excluding the boundaries
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                for (int k = 1; k <= N; k++) {
                    // Calculate the updated value for u_1[i][j][k]
                    double u_0_1_0_0 = u_1[i - 1][j][k];
                    double u_0_0_1_0 = u_1[i][j - 1][k];
                    double u_0_0_0_1 = u_1[i][j][k - 1];

                    double u_1_1_0_0 = u_1[i + 1][j][k];
                    double u_1_0_1_0 = u_1[i][j + 1][k];
                    double u_1_0_0_1 = u_1[i][j][k + 1];

                    u_1[i][j][k] = (u_1_1_0_0 + u_1_0_1_0 + u_1_0_0_1 + 
                                    u_0_1_0_0 + u_0_0_1_0 + u_0_0_0_1 + 
                                    delta * delta * f[i][j][k]) / 6.0;
                }
            }
        }
    }

    return;
}

void par_gauss_seidel(double ***u, double ***f, int N, int P)
{
    // Grid spacing
    double delta = 2.0 / N;

    // Perform P iterations
    for (int p = 0; p < P; p++) {
        #pragma omp for ordered(2) private(j, k)
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                #pragma omp ordered depend(sink: i-1, j-1) depend(sink: i-1, j) \
                                        depend(sink: i-1, j+1) depend(sink: i, j-1)
                for (int k = 1; k < N - 1; ++k) {
                    double tmp1 = (
                        u[i-1][j-1][k-1] + u[i-1][j-1][k] + u[i-1][j-1][k+1] +
                        u[i-1][j][k-1]   + u[i-1][j][k]   + u[i-1][j][k+1]   +
                        u[i-1][j+1][k-1] + u[i-1][j+1][k] + u[i-1][j+1][k+1]
                    );

                    double tmp2 = (
                        u[i][j-1][k-1] + u[i][j-1][k] + u[i][j-1][k+1] +
                        u[i][j][k-1]   + u[i][j][k]   + u[i][j][k+1]   +
                        u[i][j+1][k-1] + u[i][j+1][k] + u[i][j+1][k+1]
                    );

                    double tmp3 = (
                        u[i+1][j-1][k-1] + u[i+1][j-1][k] + u[i+1][j-1][k+1] +
                        u[i+1][j][k-1]   + u[i+1][j][k]   + u[i+1][j][k+1]   +
                        u[i+1][j+1][k-1] + u[i+1][j+1][k] + u[i+1][j+1][k+1]
                    );

                    u[i][j][k] = (tmp1 + tmp2 + tmp3) / 27.0;
                }
                
                #pragma omp ordered depend(source)
            }
        }
    }

    return;
}

// double*** gauss_seidel_omp(double ***u, double ***f, int N) {
//     double delta = 2.0 / N;
//     double delta_2 = delta * delta;
//     double factor = 1.0 / 6.0;

//     // Parallelize using OpenMP with doacross
//     #pragma omp parallel default(none) shared(u, f, N, delta_2, factor) private(i, j, k)
//     {
//         for (int i = 1; i < N + 1; ++i) {
//             #pragma omp for ordered(2) collapse(2) schedule(static, 1)
//             for (int j = 1; j < N + 1; ++j) {
//                 #pragma omp ordered depend(sink: i - 1, j - 1) \
//                                       depend(sink: i - 1, j) \
//                                       depend(sink: i - 1, j + 1) \
//                                       depend(sink: i, j - 1)
//                 for (int k = 1; k < N + 1; ++k) {
//                     u[i][j][k] = factor * (
//                         u[i - 1][j][k] + u[i + 1][j][k] +
//                         u[i][j - 1][k] + u[i][j + 1][k] +
//                         u[i][j][k - 1] + u[i][j][k + 1] +
//                         delta_2 * f[i][j][k]
//                     );
//                 }
//                 #pragma omp ordered depend(source)
//             }
//         }
//     }
//     return u;
// }