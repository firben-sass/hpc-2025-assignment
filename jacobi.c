/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>


int seq_jacobi(double *** u_0, double *** u_1, double *** f, int N, int P)
{
    double factor = 1.0/6.0;
    double delta = 2.0/N;

    for (int p = 0; p < P; p++)
    {
        for (int i = 1; i < N+1; i++)
        {
            for (int j = 1; j < N+1; j++)
            {
                for (int k = 1; k < N+1; k++)
                {
                    u_1[i][j][k] = (u_0[i-1][j][k] + u_0[i+1][j][k] + u_0[i][j-1][k] + u_0[i][j+1][k] + u_0[i][j][k-1] + u_0[i][j][k+1] + delta*delta*f[i][j][k]) * factor;
                }
            }
        }

        // Swap the pointers to update u_0 for the next iteration
        double ***temp = u_0;
        u_0 = u_1;
        u_1 = temp;
    }

    return 0;
}

int par_jacobi(double *** u_0, double *** u_1, double *** f, int N, int P)
{
    double factor = 1.0 / 6.0;
    double delta = 2.0 / N;

    for (int p = 0; p < P; p++)
    {
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 1; i < N + 1; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                for (int k = 1; k < N + 1; k++)
                {
                    u_1[i][j][k] = (u_0[i-1][j][k] + u_0[i+1][j][k] + 
                                    u_0[i][j-1][k] + u_0[i][j+1][k] + 
                                    u_0[i][j][k-1] + u_0[i][j][k+1] + 
                                    delta * delta * f[i][j][k]) * factor;
                }
            }
        }

        // Swap the pointers to update u_0 for the next iteration
        double ***temp = u_0;
        u_0 = u_1;
        u_1 = temp;
    }

    return 0;
}
