/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>
#include <stdio.h>


double*** gauss_seidel(double ***u,double ***f,int N) {
    double delta = 2.0 / (N);
    double delta_2 = delta * delta;
    double factor = 1.0 / 6.0;
    for(int i = 1; i < N + 1; i++){
        for(int j = 1; j < N + 1; j++){
           for(int k = 1; k < N + 1; k++){
                u[i][j][k] = factor*(
                u[i-1][j][k] + u[i+1][j][k] + 
                u[i][j-1][k] + u[i][j+1][k] + 
                u[i][j][k-1] + u[i][j][k+1] + 
                delta_2*f[i][j][k]
                );
            }
        }
    }
    return u;    
}

