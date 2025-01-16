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

// TODO: Spørg hjælpelærer om den er korrekt
void seq_gauss_seidel(double ***u_0, double ***u_1, double ***f, int N, int P)
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
                    double u_1_1_0_0 = u_1[i - 1][j][k];
                    double u_1_0_1_0 = u_1[i][j - 1][k];
                    double u_1_0_0_1 = u_1[i][j][k - 1];

                    double u_0_1_0_0 = u_0[i + 1][j][k];
                    double u_0_0_1_0 = u_0[i][j + 1][k];
                    double u_0_0_0_1 = u_0[i][j][k + 1];

                    u_1[i][j][k] = (u_1_1_0_0 + u_1_0_1_0 + u_1_0_0_1 + 
                                    u_0_1_0_0 + u_0_0_1_0 + u_0_0_0_1 + 
                                    delta * delta * f[i][j][k]) / 6.0;
                }
            }
        }

        // Swap the pointers to update u_0 for the next iteration
        double ***temp = u_0;
        u_0 = u_1;
        u_1 = temp;
    }

    return;
}

void par_gauss_seidel(double ***u_0, double ***u_1, double ***f, int N, int P)
{
    // Grid spacing
    double delta = 2.0 / N;

    // Perform P iterations
    for (int p = 0; p < P; p++) {
        // Parallelize the loops using OpenMP
        #pragma omp parallel for collapse(3) schedule(dynamic)
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                for (int k = 1; k <= N; k++) {
                    // Calculate the updated value for u_1[i][j][k]
                    double u_1_1_0_0 = u_1[i - 1][j][k];
                    double u_1_0_1_0 = u_1[i][j - 1][k];
                    double u_1_0_0_1 = u_1[i][j][k - 1];

                    double u_0_1_0_0 = u_0[i + 1][j][k];
                    double u_0_0_1_0 = u_0[i][j + 1][k];
                    double u_0_0_0_1 = u_0[i][j][k + 1];

                    u_1[i][j][k] = (u_1_1_0_0 + u_1_0_1_0 + u_1_0_0_1 + 
                                    u_0_1_0_0 + u_0_0_1_0 + u_0_0_0_1 + 
                                    delta * delta * f[i][j][k]) / 6.0;
                }
            }
        }

        // Swap the pointers to update u_0 for the next iteration
        double ***temp = u_0;
        u_0 = u_1;
        u_1 = temp;
    }

    return;
}
