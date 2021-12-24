/* prob.c */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "linalg.h"

double runiform() { return (double) rand() / RAND_MAX; }

double rnormal(double mean, double std) { /* Box-Muller. */
    static double cached = 0.0;
    double x, y, r, res;

    if (cached == 0.0) {
        do {
            x = 2.0 * rand() / RAND_MAX - 1;
            y = 2.0 * rand() / RAND_MAX - 1;
            r = x * x + y * y;
        } while (r == 0.0 || r > 1.0);

        double d = sqrt(-2.0 * log(r) / r);

        double n1 = x * d;
        double n2 = y * d;

        res = n1 * std + mean;
        cached = n2;
    }
    else {
        res = cached * std + mean;
        cached = 0.0;
    }

    return res;
}

double rgamma(double alpha, double beta) { /* Marsaglia and Tsang. */
    int small = FALSE;
	   double d, c, x, v, u;

    if (alpha <= 0.0 || beta <= 0.0) error("rgamma: parameters must be positive.");

    if (alpha < 1.0) { alpha += 1.0; small = TRUE; }

	   d = alpha - 1.0 / 3.0; c = 1.0 / sqrt(9.0 * d);

	   do {
		      do { x = rnormal(0.0, 1.0); v = 1.0 + c * x; } while (v <= 0.0);

	       v = CUBE(v);
	       u = runiform();
    }
    while (u > 1.0 - 0.0331 * SQUARE(SQUARE(x))
	            &&
	          log(u) > 0.5 * SQUARE(x) + d * (1.0 - v + log(v)));

    if (small) {
		      do { u = runiform(); } while (u == 0.0);
		      return pow(u, 1.0 / (alpha - 1.0)) * (d * v) / beta;
    }
    else return (d * v) / beta;
}

void rdirichlet(int n, double *a, double *x) {
    int i;
    double sum = 0.0;
    double *g = vector(n);
    for (i = 0; i < n; i++) { g[i] = rgamma(a[i], 1.0); sum += g[i]; }
    for (i = 0; i < n; i++) x[i] = (g[i] + 0.01) / sum;
    free(g);
}

void rriemann(double *cur, double *can, int n, double a0, double delta) {
	   int i;
	   double *a = vector(n);
	   for (i = 0; i < n; i++) a[i] = cur[i] * a0 * delta;
	   rdirichlet(n, a, can);
	   for (i = 0; i < n; i++) can[i] /= delta;
	   free(a);
}

void rrw(double *cur, double *can, int n, double eps, double delta) {
	   int i;
	   double tmp, sum = 0.0;
	   for (i = 0; i < n; i++) { tmp = cur[i] * exp(rnormal(0, eps) + 0.01); can[i] = tmp; sum += tmp * delta; }
	   for (i = 0; i < n; i++) can[i] = can[i] / sum;
}

double lgamma(double x) { /* NR in C. */
    int i;
    double u, v, tmp, ser;

    static const double cof[14]= {
		  57.1562356658629235      , - 59.5979603554754912      ,
          14.1360979747417471      , -  0.491913816097620199    ,
           0.339946499848118887e-4 ,    0.465236289270485756e-4 ,
        -  0.983744753048795646e-4 ,    0.158088703224912494e-3 ,
        -  0.210264441724104883e-3 ,    0.217439618115212643e-3 ,
        -  0.164318106536763890e-3 ,    0.844182239838527433e-4 ,
        -  0.261908384015814087e-4 ,    0.368991826595316234e-5
    };

    if (x <= 0.0) error("lgamma: argument must be non-negative.");

    u = v = x;

    tmp = u + 5.24218750000000000;
    tmp = (u + 0.5) * log(tmp) - tmp;
    ser = 0.999999999999997092;
    for (i = 0; i < 14; i++) ser += cof[i] / ++v;

    return tmp + log(2.5066282746310005 * ser / u);
}
