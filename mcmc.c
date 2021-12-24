/* mcmc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lib/util.h"
#include "lib/linalg.h"
#include "lib/prob.h"
#include "mcmc.h"

long main()
{
	long i, j;
	double x[N], eps, theta;

	srand(time(NULL));

    time_t start = time(NULL);

	long k = (long) ceil((T_MAX - T_MIN) / DELTA);
    double *t = vector(k + 1);
    for (i = 0; i <= k; i++) t[i] = T_MIN + i * DELTA;

    read_data(x);

    double *c = vector(k);
	for (i = 0; i < N; i++) for (j = 0; j < k; j++) if (t[j] < x[i] && x[i] < t[j + 1]) c[j]++;

    clear_results();

    for (i = 1; i <= THETA_DRAWS; i++)
    {
        theta = rgamma(ALPHA, BETA) + THETA_C;
        printf("theta = %f\n", theta);

        double *h_est = vector(k);

        metropolis_hastings(theta, h_est, &eps, t, k, c);

        save_results(h_est, k, eps);

        free(h_est);
    }

    free(t); free(c);

    printf("Elapsed: %.2f minutes.\n", difftime(time(NULL), start) / 60.0);

    return 0;
}

void read_data(double *x) {
    long i;
    FILE *fdata = fopen("data.txt", "r"); if (fdata == NULL) error("Cannot open data.txt file.");
	for (i = 0; i < N; i++) { fscanf(fdata, "%lf", &x[i]); if (x[i] < T_MIN || x[i] > T_MAX) error("Out of range data."); }
	fclose(fdata);
}

void clear_results() {
    FILE *fh_est = fopen("h_est.txt", "w");
    if (fh_est == NULL) error("Cannot create h_est.txt file.");
    fclose(fh_est);

    FILE *feps = fopen("eps.txt", "w");
    if (feps == NULL) error("Cannot create eps.txt file.");
	fclose(feps);
}

void save_results(double *h_est, long k, double eps) {
    long i;

    FILE *fh_est = fopen("h_est.txt", "a");
    if (fh_est == NULL) error("Cannot append to h_est.txt file.");
    for (i = 0; i < k - 1; i++) fprintf(fh_est, "%f\t", h_est[i]);
    fprintf(fh_est, "%f\n", h_est[k - 1]);
    fclose(fh_est);

    FILE *feps = fopen("eps.txt", "a");
    if (feps == NULL) error("Cannot append to eps.txt file.");
    fprintf(feps, "%f\n", eps);
	fclose(feps);
}

long cmpDoubles(const void *table_entry_ptr1, const void *table_entry_ptr2) {
	double x = * (double *) table_entry_ptr1;
	double y = * (double *) table_entry_ptr2;
	return (x == y) ? 0 : ((x < y) ? -1 : 1);
}

double credible_set(double *h_est, long k, double **indep, long num_hops) {
	long i, m, cr_idx;
	long j;

	double *eps = vector(num_hops);

	for (i = 0; i < num_hops; i++)
	    for (j = 0; j < k; j++)
			if (fabs(h_est[j] - indep[i][j]) > eps[i]) eps[i] = fabs(h_est[j] - indep[i][j]);

    qsort(eps, num_hops, sizeof(double), cmpDoubles);

    cr_idx = floor(CR_LEVEL * num_hops);
    double epsilon = eps[cr_idx];

	free(eps);

	return epsilon;
}

void save_chain(double *chain, long size) {
	long i;

	FILE *fchain = fopen("chain.txt", "w");
	if (fchain == NULL) error("Cannot write to chain.txt file.");
	for (i = 0; i < size; i++) fprintf(fchain, "%f\n", chain[i]);
	fclose(fchain);
}

void metropolis_hastings(double theta, double *h_est, double *ptr_eps, double *t, long k, double *c) {
    long i, j;
    long step, accepted = 0;

    double *ms = vector(k);
    double **S = matrix(k, k);
    double **L = lt_matrix(k);

    double *chain = vector(H_BURN_IN + H_DRAWS);
    long idx_mon = (long) floor(k / 2);

    long num_hops = (long) floor(H_DRAWS / H_HOP);
    double **indep = matrix(num_hops, k);
    long hop_idx = 0;

    double *h_cur = vector(k);
    for (i = 0; i < k; i++) h_cur[i] = 1.0 / (T_MAX - T_MIN);
    double *h_can = vector(k);

    for (i = 0; i < k; i++) {
		ms[i] = MEAN((t[i] + t[i + 1]) / 2);
		for (j = 0; j < k; j++) {
			S[i][j] = COV(RHO, theta, (t[i] + t[i + 1]) / 2, (t[j] + t[j + 1]) / 2);
			if (POSTERIOR) ms[i] += S[i][j] * c[j];
	    }
    }
    cholesky(S, L, k);

    for (step = 1; step < H_BURN_IN + H_DRAWS; step++)
    {
		/* rriemann(h_cur, h_can, k, A0, DELTA); */
		rrw(h_cur, h_can, k, 1e-2, DELTA);

		if (log(runiform()) < min(log_ratio(h_cur, h_can, ms, L, k), 0.0)) {
		    for (i = 0; i < k; i++) h_cur[i] = h_can[i];
		    if (step >= H_BURN_IN) accepted++;
	    }

        chain[step] = h_cur[idx_mon];

        if (step >= H_BURN_IN) {
			for (i = 0; i < k; i++) h_est[i] += (h_cur[i] / H_DRAWS);

            if (step % H_HOP == 0) {
				for (i = 0; i < k; i++) indep[hop_idx][i] = h_cur[i];
				hop_idx++;
		    }
	    }
    }

    *ptr_eps = credible_set(h_est, k, indep, num_hops);

    save_chain(chain, H_BURN_IN + H_DRAWS);

    printf("Acceptance rate = %f\n", (double) accepted / H_DRAWS);

    free(ms); free_matrix(S, k); free_lt_matrix(L, k); free(chain); free(h_cur); free(h_can); free_matrix(indep, num_hops);
}

double log_ratio(double *h_cur, double *h_can, double *ms, double **L, long k) {
    long i;
    double log_ratio = 0.0;

    for (i = 0; i < k; i++) log_ratio += ( log(h_cur[i]) - log(h_can[i]) );

    log_ratio += ( Q(h_can, ms, L, k) - Q(h_cur, ms, L, k) );

    /* log_ratio += ( log_riemann(h_cur, h_can, k) - log_riemann(h_can, h_cur, k) ); */

    return log_ratio;
}

double Q(double *h, double *ms, double **L, long k) {
    long i;
    double norm_sqr = 0.0;

    double *z = vector(k);
    for (i = 0; i < k; i++) z[i] = log(h[i]) - ms[i];

    double *u = vector(k);
    fwd_subst(L, k, u, z); /* L u = z */

    for (i = 0; i < k; i++) norm_sqr += SQUARE(u[i]);

    free(z); free(u);

    return -0.5 * norm_sqr;
}

double log_riemann(double *x, double *y, long k) {
    long i;
    double log_density = 0.0;
    double a_i;

    for (i = 0; i < k; i++) {
        a_i = A0 * DELTA * y[i];
        if (a_i > 0) log_density += a_i * log(DELTA) + (a_i - 1) * log(x[i]) - lgamma(a_i);
    }

    return log_density + lgamma(A0);
}
