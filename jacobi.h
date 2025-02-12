/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBI_H
#define _JACOBI_H

// Define prototype
double *** jacobi(double ***, double ***, double ***, int);

double *** jacobi_par(double ***, double ***, double ***, int);
#endif
