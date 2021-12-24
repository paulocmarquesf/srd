/* linalg.c */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalg.h"
#include "util.h"

double *vector(int dim) {
    int i;
    double *v = malloc(dim * sizeof(double));
    for (i = 0; i < dim; i++) v[i] = 0.0;
    return v;
}

double **matrix(int nrows, int ncols) {
    int i, j;
    double **A = malloc(nrows * ncols * sizeof(double));
    for (i = 0; i < nrows; i++) {
        A[i] = malloc(ncols * sizeof(double));
        for (j = 0; j < ncols; j++) A[i][j] = 0.0;
    }
    return A;
}

double **symm_matrix(int nrows) {
    int i, j;
    double **A = malloc((nrows * (nrows + 1) / 2) * sizeof(double));
    for (i = 0; i < nrows; i++) {
        A[i] = malloc((i + 1) * sizeof(double));
        for (j = 0; j <= i; j++) A[i][j] = 0.0;
    }
    return A;
}

double **lt_matrix(int nrows) { return symm_matrix(nrows); }

void free_symm_matrix(double **A, int nrows) {
    int i;
    for (i = 0; i < nrows; i++) free(A[i]);
    free(A);
}

void free_matrix(double **A, int nrows) { free_symm_matrix(A, nrows); }

void free_lt_matrix(double **L, int nrows) { free_symm_matrix(L, nrows); }

void print_matrix(double **A, int nrows, int ncols) {
    int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) printf("%.2e ", A[i][j]);
        printf("\n");
    }
}

void print_symm_matrix(double **A, int nrows) {
    int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < nrows; j++) {
            if (j <= i) printf("%.2e ", A[i][j]);
            else printf("%.2e ", A[j][i]);
        }
        printf("\n");
    }
}

void print_lt_matrix(double **A, int nrows) {
    int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < nrows; j++) {
            if (j <= i) printf("%.2e ", A[i][j]);
            else printf("%.2e ", 0.0);
        }
        printf("\n");
    }
}

void print_vector(double *x, int dim) {
    int i;
    for (i = 0; i < dim; i++) printf("%f ", x[i]);
    printf("\n");
}

void cholesky(double **A, double **L, int nrows) { /* A = L L' */
    int i, j, k;
    double p;

    for (i = 0; i < nrows; i++) {
        for (j = 0; j <= i; j++) {
            p = A[i][j];
            for (k = 0; k <= j - 1; k++) p -= L[i][k] * L[j][k];
            if (i == j) {
                if (p <= 0.0) error("cholesky: the matrix is not positive definite.");
                L[i][i] = sqrt(p);
            }
            else L[i][j] = p / L[j][j];
        }
    }
}

void fwd_subst(double **L, int nrows, double *x, double *b) { /* L x = b */
    int i, j;

    for (i = 0; i < nrows; i++) {
        x[i] = b[i];
        for (j = 0; j <= i - 1; j++) x[i] -= L[i][j] * x[j];
        x[i] /= L[i][i];
    }
}
