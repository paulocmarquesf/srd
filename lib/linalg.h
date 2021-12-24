/* linalg.h */

double *vector(long dim);

void print_vector(double *x, long dim);

double **matrix(long nrows, long ncols);

void print_matrix(double **A, long nrows, long ncols);

void free_matrix(double **A, long nrows);

double **symm_matrix(long nrows);

void print_symm_matrix(double **A, long nrows);

void free_symm_matrix(double **A, long nrows);

double **lt_matrix(long nrows);

void print_lt_matrix(double **A, long nrows);

void free_lt_matrix(double **L, long nrows);

void cholesky(double **A, double **L, long nrows);

void fwd_subst(double **L, long nrows, double *x, double *b);
