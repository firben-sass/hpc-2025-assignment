#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double*** gauss_seidel_omp(double ***u, double ***f, int N) {
    double delta = 2.0 / (N-1);
    double delta_2 = delta * delta;
    double factor = 1.0 / 6.0;

    // Parallelize using OpenMP with doacross
    #pragma omp parallel default(none) shared(u, f, N, delta_2, factor) private(i, j, k)
    {
        for (int i = 1; i < N + 1; ++i) {
            #pragma omp for ordered(2) collapse(2) schedule(static, 1)
            for (int j = 1; j < N + 1; ++j) {
                #pragma omp ordered depend(sink: i - 1, j - 1) \
                                      depend(sink: i - 1, j) \
                                      depend(sink: i - 1, j + 1) \
                                      depend(sink: i, j - 1)
                for (int k = 1; k < N + 1; ++k) {
                    u[i][j][k] = factor * (
                        u[i - 1][j][k] + u[i + 1][j][k] +
                        u[i][j - 1][k] + u[i][j + 1][k] +
                        u[i][j][k - 1] + u[i][j][k + 1] +
                        delta_2 * f[i][j][k]
                    );
                }
                #pragma omp ordered depend(source)
            }
        }
    }
    return u;
}
