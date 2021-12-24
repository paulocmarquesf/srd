/* util.h */

#include <time.h>

#define TRUE  1
#define FALSE 0

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

#define SQUARE(x) ((x) * (x))
#define CUBE(x)   ((x) * (x) * (x))

void error(char *msg);

void print_status(long step, time_t start, long h_burn_in, long h_draws);
