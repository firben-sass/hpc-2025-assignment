#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// int gauss_seidel_omp(double ***u,double ***f,int N, int max_iter, double tolerance) {
//     double u_old;
//     double delta = 2.0 / (N-1);
//     double delta_2 = delta * delta;
//     double factor = 1.0 / 6.0;
//     double diff = tolerance + 1.0;
//     int iter = 0;
//     while(iter < max_iter && diff > tolerance){
//         diff = 0.0;
//         #pragma omp for ordered(2) private(j,k)
//         for(int i = 1; i < N + 1; i++){
//             for(int j = 1; j < N + 1; j++){
//                  #pragma omp ordered depend(sink: i-1,j-1) depend(sink: i-1,j) \
//                 depend(sink: i-1,j+1) depend(sink: i,j-1)
//                 for(int k = 1; k < N + 1; k++){
//                     u_old = u[i][j][k];
//                     u[i][j][k] = factor*(
//                         u[i-1][j][k] + u[i+1][j][k] + 
//                         u[i][j-1][k] + u[i][j+1][k] + 
//                         u[i][j][k-1] + u[i][j][k+1] + 
//                         delta_2*f[i][j][k]
//                     );
//                     diff += (u[i][j][k] - u_old) * (u[i][j][k] - u_old);
//                 }
//                 #pragma omp ordered depend(source)
//             }
//         }
//         iter++;
//     }
//     return(iter);
// }


// int gauss_seidel_omp(double ***u,double ***f,int N, int max_iter, double tolerance) {
//     double u_old;
//     double delta = 2.0 / (N-1);
//     double delta_2 = delta * delta;
//     double factor = 1.0 / 6.0;
//     double diff = tolerance + 1.0;
//     int iter = 0;
//     while(iter < max_iter && diff > tolerance){
//     diff = 0.0;  // Thread-local variable for reduction
    
//     #pragma omp parallel
//         {
//             // We break the iteration into smaller chunks
//             #pragma omp for schedule(static, 1)
//             for (int i = 1; i < N - 1; i++) {
//                 for (int j = 1; j < N - 1; j++) {
//                     for (int k = 1; k < N - 1; k++) {
//                         u_old = u[i][j][k];
                        
//                         // Apply the Gauss-Seidel update formula
//                         u[i][j][k] = factor * (
//                             u[i-1][j][k] + u[i+1][j][k] +
//                             u[i][j-1][k] + u[i][j+1][k] +
//                             u[i][j][k-1] + u[i][j][k+1] +
//                             delta_2 * f[i][j][k]
//                         );

//                         // Accumulate the difference
//                         diff += (u[i][j][k] - u_old) * (u[i][j][k] - u_old);
//                     }
//                 }
//             }
//         }
    
        
//         iter++;
//     }
//     return(iter);
// }


// int gauss_seidel_omp(double ***u, double ***f, int N, int max_iter, double tolerance) {
//     double u_old;
//     double delta = 2.0 / (N - 1);
//     double delta_2 = delta * delta;
//     double factor = 1.0 / 6.0;
//     double diff = tolerance + 1.0;
//     int iter = 0;

//     while (iter < max_iter && diff > tolerance) {
//         diff = 0.0;

//         #pragma omp parallel
//         {
//             #pragma omp for ordered(2) private(j, k, u_old) reduction(+:diff)
//             for (int i = 1; i < N + 1; i++) {
//                 for (int j = 1; j < N + 1; j++) {
//                     #pragma omp ordered depend(sink: i - 1, j - 1) depend(sink: i - 1, j) \
//                                         depend(sink: i - 1, j + 1) depend(sink: i, j - 1)
//                     for (int k = 1; k < N +1; k++) {
//                         u_old = u[i][j][k];
//                         u[i][j][k] = factor * (
//                             u[i - 1][j][k] + u[i + 1][j][k] +
//                             u[i][j - 1][k] + u[i][j + 1][k] +
//                             u[i][j][k - 1] + u[i][j][k + 1] +
//                             delta_2 * f[i][j][k]
//                         );
//                         diff += (u[i][j][k] - u_old) * (u[i][j][k] - u_old);
//                     }
//                     #pragma omp ordered depend(source)
//                 }
//             }
//         }
//         diff = sqrt(diff);
//         iter++;
//     }

//     return iter;
// }


// int gauss_seidel_omp(double ***u, double ***f, int N,int max_iter, double tolerance) {
//     double delta = 2.0 / (N - 1);
//     double delta_2 = delta * delta;
//     double factor = 1.0 / 6.0;
//     int iter;
//     for (iter = 0; iter < max_iter; iter++) {
//         #pragma omp for ordered(2)private(j,k)
//         for (int i = 1; i < N - 1; i++) {
//             for (int j = 1; j < N - 1; j++) {
//                 #pragma omp ordered depend(sink: i-1,j-1) depend(sink: i-1,j) \
//                                     depend(sink: i-1,j+1) depend(sink: i,j-1)
//                 for (int k = 1; k < N - 1; k++) {
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
//     return iter;
// }

// double*** gauss_seidel_omp(double ***u, double ***f, int N) {
//     double delta = 2.0 / N;
//     double delta_2 = delta * delta;
//     double factor = 1.0 / 6.0;

//     // Parallelize using OpenMP with doacross
//      #pragma omp parallel for collapse(3)
//     for (int i = 1; i < N + 1; i++) {
//     for (int j = 1; j < N + 1; j++) {
//         for (int k = 1; k < N + 1; k++) {
//             if ((i + j + k) % 2 == 0) { // Red update
//                 u[i][j][k] = factor * (
//                     u[i - 1][j][k] + u[i + 1][j][k] +
//                     u[i][j - 1][k] + u[i][j + 1][k] +
//                     u[i][j][k - 1] + u[i][j][k + 1] +
//                     delta_2 * f[i][j][k]
//                 );
//             }
//         }
//     }
// }

//     #pragma omp parallel for collapse(3)
//     for (int i = 1; i < N + 1; i++) {
//         for (int j = 1; j < N + 1; j++) {
//             for (int k = 1; k < N + 1; k++) {
//                 if ((i + j + k) % 2 == 1) { // Black update
//                     u[i][j][k] = factor * (
//                         u[i - 1][j][k] + u[i + 1][j][k] +
//                         u[i][j - 1][k] + u[i][j + 1][k] +
//                         u[i][j][k - 1] + u[i][j][k + 1] +
//                         delta_2 * f[i][j][k]
//                     );
//                 }
//             }
//         }
//     }
//     return u;
// }

double*** gauss_seidel_omp(double ***u, double ***f, int N) {
    double delta = 2.0 / N;
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