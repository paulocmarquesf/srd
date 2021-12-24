/* linalg.c */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalg.h"

double *vector(long dim) {
    long i;
 	double *v = malloc(dim * sizeof(double));
	for (i = 0; i < dim; i++) v[i] = 0.0;
	return v;
}

double **matrix(long nrows, long ncols) {
	long i, j;
    double **A = malloc(nrows * ncols * sizeof(double));
    for (i = 0; i < nrows; i++) {
		A[i] = malloc(ncols * sizeof(double));
		for (j = 0; j < ncols; j++) A[i][j] = 0.0;
    }
    return A;
}

double **symm_matrix(long nrows) {
	long i, j;
    double **A = malloc((nrows * (nrows + 1) / 2) * sizeof(double));
    for (i = 0; i < nrows; i++) {
		A[i] = malloc((i + 1) * sizeof(double));
		for (j = 0; j <= i; j++) A[i][j] = 0.0;
    }
    return A;
}

double **lt_matrix(long nrows) { return symm_matrix(nrows); }

void free_symm_matrix(double **A, long nrows) {
    long i;
	for (i = 0; i < nrows; i++) free(A[i]);
	free(A);
}

void free_matrix(double **A, long nrows) { free_symm_matrix(A, nrows); }

void free_lt_matrix(double **L, long nrows) { free_symm_matrix(L, nrows); }

void print_matrix(double **A, long nrows, long ncols) {
	long i, j;
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) printf("%.2e ", A[i][j]);
		printf("\n");
    }
}

void print_symm_matrix(double **A, long nrows) {
	long i, j;
	for (i = 0; i < nrows; i++) {
	    for (j = 0; j < nrows; j++) {
			if (j <= i) printf("%.2e ", A[i][j]);
			else printf("%.2e ", A[j][i]);
	    }
	    printf("\n");
    }
}

void print_lt_matrix(double **A, long nrows) {
	long i, j;
	for (i = 0; i < nrows; i++) {
	    for (j = 0; j < nrows; j++) {
			if (j <= i) printf("%.2e ", A[i][j]);
			else printf("%.2e ", 0.0);
	    }
	    printf("\n");
    }
}

void print_vector(double *x, long dim) {
	long i;
	for (i = 0; i < dim; i++) printf("%f ", x[i]);
	printf("\n");
}

void cholesky(double **A, double **L, long nrows) { /* A = L L' */
    long i, j, k;
    double p;

    for (i = 0; i < nrows; i++) {
		for (j = 0; j <= i; j++) {
			p = A[i][j];
			for (k = 0; k <= j - 1; k++) p -= L[i][k] * L[j][k];
			if (i == j) {
		        if (p <= 0.0)
		            error("cholesky: the matrix is not positive definite.");
				L[i][i] = sqrt(p);
		    }
		    else L[i][j] = p / L[j][j];
	    }
    }
}

void fwd_subst(double **L, long nrows, double *x, double *b) { /* L x = b */
    long i, j;

    for (i = 0; i < nrows; i++) {
        x[i] = b[i];
        for (j = 0; j <= i - 1; j++) x[i] -= L[i][j] * x[j];
        x[i] /= L[i][i];
    }
}
