/*****************************************************************************
 * This file is a part of gEcon.                                             *
 *                                                                           *
 * (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015 *
 * (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018   *
 * License terms can be found in the file 'LICENCE'                          *
 *                                                                           *
 * Author: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         *
 *****************************************************************************/

/** \file zgges_interface.c
 * \brief Function zgges_interface computes qz decomposition
 * using Lapack ZGGES routine
 */

#include <math.h>
#include <stdio.h>
#include <complex.h>

extern void zgges_(const char *jobvsl,
            const char *jobvsr,
            const char *sort,
            const int *selectg,
            int *n1,
            double complex *A,
            int *n2,
            double complex  *B,
            int *n3,
            const int *sdim,
            double complex  *ALPHA,
            double complex  *BETA,
            double complex  *VSL,
            int *n4,
            double complex *VSR,
            int *n5,
            double complex (*WORK)[],
            int *lwork,
            double (*RWORK)[],
            int (*BWORK)[],
            int **info);

void
zgges_interface(int *size,
                int *info,
                double complex *A,
                double complex *B,
                double complex *ALPHA,
                double complex *BETA,
                double complex *VSL,
                double complex *VSR)
{
    // Lapack routine settings
    const char jobvsl = 'V';
    const char jobvsr = 'V';
    const char sort = 'N';
    const int selectg = 0;
    int n = *size;
    const int sdim = 0;

    // Allocating workspace
    int lwork = 3 * *size;
    double complex WORK[lwork];
    double RWORK[lwork] ;
    int BWORK[lwork];

    // Calling Fortran routine zgges_
    zgges_( &jobvsl,  /* (input) CHARACTER*1
                      'N':  do not compute the left Schur vectors;
                      'V': compute the left Schur vectors. */
            &jobvsr,  /* (input) CHARACTER*1
                      'N':  do not compute the right Schur vectors;
                      'V':  compute the right Schur vectors. */
            &sort,    /* (input) CHARACTER*1
                      'N':  Eigenvalues are not ordered;
                      'S':  Eigenvalues are ordered (see SELCTG)  */
            &selectg, /* (external procedure) LOGICAL FUNCTION
                       An eigenvalue ALPHA(j)/BETA(j) is selected
                       if SELCTG(ALPHA(j),BETA(j)) is true. */
            &n,       /*  (input) INTEGER The order of the matrices
                       A, B, VSL, and VSR.  N >= 0. */
            A,        /* (input/output) COMPLEX*16 array,
                       dimension (LDA, N) On exit,
                       A has been overwritten by its
                       generalized Schur S */
            &n,       /* (input) INTEGER  The leading dimension of A.
                       LDA >= max(1,N). */
            B,        /*  (input/output) COMPLEX*16 array,
                       dimension (LDA, N)  On exit,
                       B has been overwritten
                       by its generalized Schur T */
            &n,       /*  (input) INTEGER  The leading dimension of B.
                       LDB >= max(1,N). */
            &sdim,    /*  (output) INTEGER
                       If SORT = 'S', SDIM = number of eigenvalues
                       (after sorting) */
            ALPHA,    /* (output) COMPLEX*16 array,
                       dimension (N) generalized eigenvalues.
                       ALPHA(j), j=1,...,N  and  BETA(j), j=1,...,N
                       are the diagonals of the complex Schur form (A,B) */
            BETA,
            VSL,      /* (output) COMPLEX*16 array,
                       dimension (LDVSL,N) If JOBVSL = 'V',
                       VSL will contain the left Schur vectors. */
            &n,       /* (input) INTEGER
                       The leading dimension of the matrix VSL. */
            VSR,      /* (output) COMPLEX*16 array, dimension (LDVSL,N)
                       If JOBVSR = 'V',
                       VSR will contain the right Schur vectors. */
            &n,       /* (input) INTEGER The leading dimension
                       of the matrix VSR. */
            &WORK,    /*    (workspace/output) COMPLEX*16 array,
                       dimension (MAX(1,LWORK)) */
            &lwork,   /* The dimension of the array WORK.
                       LWORK >= max(1,2*N). */
            &RWORK,   /* (workspace) DOUBLE PRECISION array,
                       dimension (N) */
            &BWORK,   /* (workspace) LOGICAL array, dimension (N) */
            &info );  /*   info integer */
}
