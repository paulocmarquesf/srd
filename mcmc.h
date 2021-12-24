/* mcmc.h */

#define POSTERIOR TRUE

#define CR_LEVEL 0.95

#define N 900000 /* sample size */

#define T_MIN -0.5
#define T_MAX  0.5
#define DELTA  0.01

#define MEAN(x) 1.0

#define COV(rho, theta, x, y) ((rho) * exp(-(theta) * SQUARE((x) - (y))))

#define THETA_C 20000.0
#define ALPHA 2.0
#define BETA 1.0
#define THETA_DRAWS 2

#define RHO 0.05

#define H_BURN_IN  10000
#define H_DRAWS    100000
#define H_HOP      10
#define H_MON_STEP 1

#define A0 1.0e6

void read_data(double *x);

void clear_results();

void save_results(double *h_est, long k, double eps);

long cmpDoubles(const void *table_entry_ptr1, const void *table_entry_ptr2);

double credible_set(double *h_est, long k, double **indep, long num_hops);

void save_chain(double *chain, long size);

void metropolis_hastings(double theta, double *h_est, double *ptr_eps, double *t, long k, double *c);

double log_ratio(double *h_cur, double *h_can, double *ms, double **L, long k);

double Q(double *h, double *ms, double **L, long k);

double log_riemann(double *x, double *y, long k);
