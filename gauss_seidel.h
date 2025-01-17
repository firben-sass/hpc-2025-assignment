/* gauss_seidel.h - Poisson problem
 *
 */
#ifndef _GAUSS_SEIDEL_H
#define _GAUSS_SEIDEL_H

// #pragma message("Including gauss_seidel.h")

void _seq_rec_gauss_seidel(double ***, double ***, double ***, int, int, int, int);

void seq_rec_gauss_seidel(double ***, double ***, double ***, int, int);

void seq_gauss_seidel(double ***, double ***, int, int);

void par_gauss_seidel(double ***, double ***, double ***, int, int);

// double*** gauss_seidel_omp(double ***, double ***, int)

#endif
