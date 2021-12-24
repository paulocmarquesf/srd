/* prob.h */

double runiform();

double rnormal(double mean, double std);

double rgamma(double alpha, double beta);

void rdirichlet(int n, double *a, double *x);

void rriemann(double *cur, double *can, int n, double a0, double delta);

void rrw(double *cur, double *can, int n, double eps, double delta);

double lgamma(double x);
