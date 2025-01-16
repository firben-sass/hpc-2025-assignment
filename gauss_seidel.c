/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

void seq_gauss_seidel(double *** u_0, double *** u_1, double *** f, int N, int P)
{
    for (int p; p < P; p++)
    {
        _seq_gauss_seidel(u_0, u_1, f, N, N+1, N+1, N+1);
        u_0 = u_1;
    }
    
    return;
}

void _seq_gauss_seidel(double *** u_0, double *** u_1, double *** f, int N, int i, int j, int k)
{   
    // If in a wall
    if (i == 0 || i == N+1 || j == 0 || j == N+1 || k == 0 || k == N+1)
        return;

    _seq_gauss_seidel(u_0, u_1, f, N, i-1, j, k);
    double u_1_1_0_0 = u_1[i-1][j][k];
    _seq_gauss_seidel(u_0, u_1, f, N, i, j-1, k);
    double u_1_0_1_0 = u_1[i][j-1][k];
    _seq_gauss_seidel(u_0, u_1, f, N, i, j, k-1);
    double u_1_0_0_1 = u_1[i][j][k-1];

    double u_0_1_0_0 = u_0[i+1][j][k];
    double u_0_0_1_0 = u_0[i][j+1][k];
    double u_0_0_0_1 = u_0[i][j][k+1];

    double delta = 2.0/N;
    u_1[i][j][k] = (u_1_1_0_0 + u_1_0_1_0 + u_1_0_0_1 + u_0_1_0_0 + u_0_0_1_0 + u_0_0_0_1 + delta*delta*f[i][j][k]) / 6.0;

    return;
}