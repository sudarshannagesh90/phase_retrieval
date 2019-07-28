#include "mex.h"
#include <string.h>
#include <math.h>
#include <complex.h>
#include "gsl/gsl_sf_bessel.h"
#include <time.h>
#include <stdlib.h>

/* SOLVERS */
void cgamp (
        size_t n, size_t m, double *y, int *ir, int *jc,
        double delta, double *pr_prmts,
        int t_max, double eps, double damp, double vnf, int disp, FILE *output,
        complex double *a, double *c
    );

/* PRIORS */
void prior( complex double r, double sig, double *p_prms, complex double *a, double *c, double vnf);

/* CHANNELS */
void channel( double y, complex double o, double v, double delta,
        complex double *g, double *dg, double vnf );


/* COMMON */
void sort_rand( int n, int *seq );
double drand(double dMin, double dMax);

static inline double max( double a, double b ) { return a > b ? a : b; }
static inline double min( double a, double b ) { return a < b ? a : b; }
static inline int rand_int( int k ) { return rand() / (RAND_MAX / k + 1); }
