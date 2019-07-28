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

#include "cswgamp.h"

/* Interface b/w MATLAB and C */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *ry;
    double *p_prms, *rinit, *cinit, t_max, eps, damp, vnf;
    double *y,*ra, *ca, *c, delta, *diff_y, *dy;
    complex double *a;

    FILE *output;

    int *ir, *jc;

    mwIndex *ir_, *jc_;
    mxArray *opt_tm, *opt_ep,*opt_del, *opt_prms, *opt_in, *opt_da,
            *opt_di, *opt_vn, *opt_ou;

    size_t m, n, nnz;
    unsigned int mu, i, key;
    int disp;

    ry = mxGetPr(prhs[0]); //cy = mxGetPi(prhs[0]);
    //rF = mxGetPr(prhs[1]); //cF = mxGetPi(prhs[1]);
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    diff_y = malloc(sizeof(double) * t_max);

    if (nrhs > 2) {
        opt_del = mxGetField(prhs[2], 0, "delta");
        opt_prms = mxGetField(prhs[2], 0, "priorPrms");
        opt_tm = mxGetField(prhs[2], 0, "maxIter");
        opt_ep = mxGetField(prhs[2], 0, "prec");
        opt_in = mxGetField(prhs[2], 0, "initState");
        opt_da = mxGetField(prhs[2], 0, "damp");
        opt_di = mxGetField(prhs[2], 0, "display");
        opt_ou = mxGetField(prhs[2], 0, "output");
        opt_vn = mxGetField(prhs[2], 0, "vnf");

    } else {
        opt_tm = opt_ep = opt_del = opt_prms = opt_in = opt_da = opt_vn = opt_di = opt_ou = NULL;
    }
                     /* Prior */
    p_prms = malloc(sizeof(complex double) * 7);
    if (opt_prms)
        p_prms = mxGetPr(opt_prms);
    else {
        p_prms[0] = 0.1;//rho
        p_prms[1] = 0.5;//gb_mean
        p_prms[2] = 1e-5;//gb_var
    }

    //for(i=n*p_prms[5]-1;i>(n-1)*p_prms[5];i--)
    //    printf("i=%d; %.4e\n", i, W[i]);

    t_max = opt_tm ? *mxGetPr(opt_tm) : 250;
    eps = opt_ep ? *mxGetPr(opt_ep) : 1e-13;
    delta = opt_del ? *mxGetPr(opt_del) : 1e-8;
    damp = opt_da ? *mxGetPr(opt_da) : 0;
    vnf = opt_vn ? *mxGetPr(opt_vn) : 0.5;
    disp = opt_di ? *mxGetPr(opt_di) : 1;
    output = opt_ou ? fopen(mxArrayToString(opt_ou), "w") : NULL;


    /* Set-up output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);

    ra = mxGetPr(plhs[0]); ca = mxGetPi(plhs[0]);
    c = mxGetPr(plhs[1]);

    if (opt_in) {
        rinit = mxGetPr(opt_in); cinit = mxGetPi(opt_in);
        for (i = 0; i < n; i++) {
            ra[i] = rinit[i], c[i] = rinit[n + i];
            ca[i] = cinit ? cinit[i] : 0.;
        }
    } else {
        for (i = 0; i < n; i++) ra[i] = 0., ca[i] = 0., c[i] = 1.;
    }


    /* Run algorithm */
    if (mxIsSparse(prhs[1])) {
        /* Generate arrays w/ indexes */
        ir_ = mxGetIr(prhs[1]);
        jc_ = mxGetJc(prhs[1]);
        nnz = jc_[n];

        ir = mxMalloc(sizeof(int) * nnz);
        jc = mxMalloc(sizeof(int) * (n + 1));
        for (key = 0; key < nnz; key++) ir[key] = ir_[key];
        for (i = 0; i < n + 1; i++) jc[i] = jc_[i];

        y = mxMalloc(sizeof(double) * m);//complex
        a = mxMalloc(sizeof(complex double) * n);

        for (mu = 0; mu < m; mu++) y[mu] = ry[mu];// + (cy ? cy[mu] * I : 0);
        //for (key = 0; key < nnz; key++) F[key] = rF[key];// + (cF ? cF[key] * I : 0);
        for (i = 0; i < n; i++) a[i] = ra[i] + (ca ? ca[i] * I : 0);

        cgamp(n, m, y, ir, jc,
            delta, p_prms,
            t_max, eps, damp, vnf, disp, output,
            a, c);



        for (i = 0; i < n; i++) ra[i] = creal(a[i]);
        if (ca)
            for (i = 0; i < n; i++) ca[i] = cimag(a[i]);



    } else {
        mexErrMsgTxt("Convert F to sparse before using this solver.");
    };

    /* Dealloc. structures */
    mxFree(a);
    mxFree(y);
    mxFree(jc);
    mxFree(ir);

    if (output) fclose(output);
}
