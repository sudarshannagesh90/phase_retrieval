/*
 * Copyright 2016 Boshra Rajaei <b.rajaei@sadjad.ac.ir>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../cswgamp.h"
//#include <time.h>

/* (Sequential) AMP for sparse matrices */
void cgamp (
        size_t n, size_t m, double *y, int *ir, int *jc,
        double delta, double *p_prms,
        int t_max, double eps, double damp, double vnf, int disp, FILE *output,
        complex double *a, double *c
    ) {
    complex double *a_proj, *o, *g, *g_old, *r;
    complex double a_old, g_proj;
    double *c_proj, *v, *dg, *sig, *rhos;
    double c_old, dg_proj, v_old;
    double diff, diff_y;
    //time_t t1, t2;

    unsigned int i, mu, idx, t;
    int *seq, key;

    /* Alloc. structures */
    a_proj = malloc(sizeof(complex double) * m);
    o = malloc(sizeof(complex double) * m);
    g = malloc(sizeof(complex double) * m);
    g_old = malloc(sizeof(complex double) * m);
    c_proj = malloc(sizeof(double) * m);
    v = malloc(sizeof(double) * m);
    dg = malloc(sizeof(double) * m);
    seq = malloc(sizeof(int) * n);
    r = malloc(sizeof(complex double) * n);
    sig = malloc(sizeof(double) * n);


    srand(time(NULL));
    
    if(disp){
        printf("delta: %.2e; t_max: %3d; eps: %.2e; damp %.2e; vnf: %.2e\n", delta, t_max, eps, damp, vnf);
        printf("c[0]: %.2e; a[0]: %.2e\n", c[0], a[0]);
        printf("prmts: %.2e; %.2e; %.2e;\n", p_prms[0], p_prms[1], p_prms[2]);
    }

    if (!a_proj || !c_proj || !o || !v || !g || !g_old || !dg || !seq)
        mexErrMsgTxt("Failure in allocating memory.");

    for (mu = 0; mu < m; mu++) g[mu] = 0;

    //t1 = time(NULL);
    for (t = 0; t < t_max; t++) {
        /* Generate random permutation */
        for (key = 0; key < n; key++) seq[key] = key;
        sort_rand(n, seq);

        /* Update a_proj and c_proj */
        for (mu = 0; mu < m; mu++)
            a_proj[mu] = c_proj[mu] = 0;
        for (i = 0; i < n; i++) {
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += a[i];
                c_proj[ ir[idx] ] += c[i];
            }
        }

        /* Update w and v */
        for (mu = 0; mu < m; mu++) {
            o[mu] = a_proj[mu] - c_proj[mu] * g[mu];
            v[mu] = c_proj[mu];
        }

        for (mu = 0; mu < m; mu++)
            channel(y[mu], o[mu], v[mu], delta, &g[mu], &dg[mu], vnf);

        for (mu = 0; mu < m; mu++) g_old[mu] = g[mu];


        /* Sweep over all n variables, in random order */
        diff = 0.;

        for (key = 0; key < n; key++) {
            i = seq[key];

            /* Update {sig, r}, {a, c} */
            a_old = a[i], c_old = c[i];

            g_proj = dg_proj = 0.; /* Dot products: F * g and -F^2 * dg */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                g_proj += g[ ir[idx] ];
                dg_proj -= dg[ ir[idx] ];
            }
            if(t>0)
                sig[i] = damp * sig[i] + (1 - damp) * (1. / dg_proj);
            else
                sig[i] = 1. / dg_proj;

            if(sig[i] < 0) sig[i] = 0.1*vnf;

            if(t>0)
                r[i] = damp * r[i] + (1 - damp) * (a[i] + sig[i] * g_proj);
            else
                r[i] = a[i] + sig[i] * g_proj;

            prior(r[i], sig[i], p_prms, &a[i], &c[i], vnf);


            /* Update {w, v}, {g, dg} */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                mu = ir[idx];
                v_old = v[mu];

                v[mu] += (c[i] - c_old);
                o[mu] += (a[i] - a_old)
                        - (v[mu] - v_old) * g_old[mu];

                channel(y[mu], o[mu], v[mu], delta, &g[mu], &dg[mu], vnf);
            }

            diff += cabs(a[i] - a_old);
        }
        diff /= n;

        for (mu = 0; mu < m; mu++)
            a_proj[mu] = 0;
        for (i = 0; i < n; i++) {
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += a[i];
            }
        }
        diff_y = 0;
        for (mu = 0; mu < m; mu++)
            diff_y += abs(y[mu] - cabs(a_proj[mu]));
        diff_y /= m;

        if (output)
            fprintf(output, "%g\n", diff_y);
        mexEvalString("drawnow");

        /* Print some info. */
        if (disp)
            printf("t: %3d; diff: %.4e; diff_y: %.4e\n", t, diff, diff_y);
        //mexEvalString("drawnow");

        /* Check for convergence */
        if (diff < eps || isnan(diff)) break;

    }

    //t2 = time(NULL);
    //printf("execution time: %d\n", t2-t1);

    /* Dealloc. structures */
    free(seq);
    free(dg);
    free(g_old);
    free(g);
    free(v);
    free(o);
    free(c_proj);
    free(a_proj);
    free(r);
    free(sig);
}
