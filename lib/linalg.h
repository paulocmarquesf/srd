/* linalg.h */

double *vector(int dim);

void print_vector(double *x, int dim);

double **matrix(int nrows, int ncols);

void print_matrix(double **A, int nrows, int ncols);

void free_matrix(double **A, int nrows);

double **symm_matrix(int nrows);

void print_symm_matrix(double **A, int nrows);

void free_symm_matrix(double **A, int nrows);

double **lt_matrix(int nrows);

void print_lt_matrix(double **A, int nrows);

void free_lt_matrix(double **L, int nrows);

void cholesky(double **A, double **L, int nrows);

void fwd_subst(double **L, int nrows, double *x, double *b);
