/* util.c */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "util.h"

void error(char *msg) { printf("\n%s\n", msg); abort(); }

void print_status(int step, time_t start, int h_burn_in, int h_draws) {
    double p = 100.0 * ((double) step / (h_burn_in + h_draws));
    double elapsed = difftime(time(NULL), start);
    printf("Iteration %i of %i (%.2f\%). Remaining %.2f minutes.\n",
           step, h_burn_in + h_draws, p, ((100 - p) / p) * elapsed / 60);
    fflush(stdout);
}
