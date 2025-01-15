/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>


int jacobi(double *** u_0, double *** u_1, double *** f, double Delta, int N, int P)
{
    for (int p = 0; p < P; p++)
    {
        for (int i = 1; i < N+1; i++)
        {
            for (int j = 1; j < N+1; j++)
            {
                for (int k = 1; k < N+1; k++)
                {
                    u_1[i][j][k] = (u_0[i-1][j][k] + u_0[i+1][j][k] + u_0[i][j-1][k] + u_0[i][j+1][k] + u_0[i][j][k-1] + u_0[i][j][k+1] + Delta*Delta*f[i][j][k]) / 6;
                }
            }
        }

        u_0 = u_1;
    }

    return 0;
}
