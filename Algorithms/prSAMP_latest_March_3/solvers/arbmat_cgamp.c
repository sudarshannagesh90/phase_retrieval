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

#include "../arbmat_cswgamp.h"
//#include <time.h>

/* (Sequential) AMP for sparse matrices */
void arbmat_cgamp (
        size_t n, size_t m, double *y, complex double *F, int *ir, int *jc,
        double delta, double *p_prms,
        int t_max, int par_frac, double eps, double damp_seq, double damp_par, double vnf, int disp, FILE *output,
        complex double *a, double *c
    ) {
    complex double *a_proj, *a_old_vec, *o, *g, *g_old, *r, *dg_proj_vec, *g_test, *a_test;
    complex double a_old, g_proj;
    double *c_proj, *v, *dg, *sig, *rhos, *dg_test, *c_test;
    double c_old, dg_proj, v_old;
    double diff, diff_y;
	double temp;
	double a_norm, r_norm, sig_ave, x_var_init, dg_ave;
	int myChar;
    //time_t t1, t2;

    unsigned int i, mu, idx, t;
    int *seq, key;

    /* Alloc. structures */
    a_test = malloc(sizeof(complex double) * n);
    c_test = malloc(sizeof(double) * n);;
    g_test = malloc(sizeof(complex double) * m);
    dg_test = malloc(sizeof(complex double) * m);
    dg_proj_vec = malloc(sizeof(complex double) * n);
	a_old_vec = malloc(sizeof(complex double) * n);
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
        printf("delta: %.2e; t_max: %3d; eps: %.2e; damp_par %.2e; vnf: %.2e\n", delta, t_max, eps, damp_par, vnf);
        printf("c[0]: %.2e; a[0]: %.2e\n", c[0], a[0]);
        printf("prmts: %.2e; %.2e; %.2e;\n", p_prms[0], p_prms[1], p_prms[2]);
    }

    if (!a_proj || !c_proj || !o || !v || !g || !g_old || !dg || !seq)
        mexErrMsgTxt("Failure in allocating memory.");

    for (mu = 0; mu < m; mu++) g[mu] = 0;
	
	a_norm=0.0;
	for (i = 0; i < n; i++){
		a_norm+=pow(cabs(a[i]),2);
		//printf("pow(cabs(a[i]),2) %.4e \n", pow(cabs(a[i]),2));
	}
	a_norm=sqrt(a_norm);
	printf("Norm of a_o is %.4e \n \n", a_norm);
    
    x_var_init=1.0/n;//oracle variance of x.  Not quite accruate because x is complex.  Will need tpo play with this later.
    /* Try initializitg sig */
    for (i = 0; i < n; i++){
//         sig[i]=x_var_init;
        sig[i]=0.0;
        r[i]=a[i];
	}
    
    
    sig_ave=0.0;
    for (i = 0; i < n; i++)
        sig_ave+=sig[i];
    sig_ave=sig_ave/n;
    printf("Average Sig is %.4e \n", sig_ave);

    r_norm=0.0;
    for (i = 0; i < n; i++)
        r_norm+=pow(cabs(r[i]),2);
    r_norm=sqrt(r_norm);
    printf("Norm of R is %.4e \n", a_norm);

    a_norm=0.0;
    for (i = 0; i < n; i++)
        a_norm+=pow(cabs(a[i]),2);
    a_norm=sqrt(a_norm);
    printf("Norm of Estimate is %.4e \n \n", a_norm);
    
    printf("a[0] is %f + i%f\n", creal(a[0]), cimag(a[0]));
    printf("c[0] is %.4e \n", c[0]);
    printf("g[0] is %.4e \n", g[0]);
    //printf("dg[1] is %.4e \n", dg[1]);
    printf("r[0] is %f + i%f\n", creal(r[0]), cimag(r[0]));
    printf("sig[0] is %.4e \n \n", sig[0]);
    

    //t1 = time(NULL);
    for (t = 0; t < t_max; t++) {
		if ((t%par_frac)!=(par_frac-1)){//Perform sequential updates on all but 1/par_frac of the updates
			/* Generate random permutation */
			for (key = 0; key < n; key++) seq[key] = key;
			sort_rand(n, seq);

			/* Update a_proj and c_proj */
			for (mu = 0; mu < m; mu++)
				a_proj[mu] = c_proj[mu] = 0;
			for (i = 0; i < n; i++) {
				for (idx = jc[i]; idx < jc[i + 1]; idx++) {
					a_proj[ ir[idx] ] += F[idx] * a[i];
					c_proj[ ir[idx] ] += (conj(F[idx]) * F[idx]) * c[i];
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
			printf("Sequential Update: \n");
			for (key = 0; key < n; key++) {
				i = seq[key];

				/* Update {sig, r}, {a, c} */
				a_old = a[i], c_old = c[i];

				g_proj = dg_proj = 0.; /* Dot products: F * g and -F^2 * dg */
				for (idx = jc[i]; idx < jc[i + 1]; idx++) {
					g_proj += F[idx] * g[ ir[idx] ];
					dg_proj -= (conj(F[idx]) * F[idx]) * dg[ ir[idx] ];
				}
				if(t>0)
					sig[i] = damp_seq * sig[i] + (1 - damp_seq) * (1. / dg_proj);
				else
					sig[i] = 1. / dg_proj;

				if(sig[i] < 0) sig[i] = 0.1*vnf;

				if(t>0)
					r[i] = damp_seq * r[i] + (1 - damp_seq) * (a[i] + sig[i] * g_proj);
				else
					r[i] = a[i] + sig[i] * g_proj;

				// prior(r[i], sig[i], p_prms, &a[i], &c[i], vnf);
				prior_nzmean_i(n, i, r[i], sig[i], p_prms, &a[i], &c[i], vnf);

				/* Update {w, v}, {g, dg} */
				for (idx = jc[i]; idx < jc[i + 1]; idx++) {
					mu = ir[idx];
					v_old = v[mu];

					v[mu] += (conj(F[idx]) * F[idx]) * (c[i] - c_old);
					o[mu] += F[idx] * (a[i] - a_old)
							- (v[mu] - v_old) * g_old[mu];

					channel(y[mu], o[mu], v[mu], delta, &g[mu], &dg[mu], vnf);
				}

				diff += cabs(a[i] - a_old);
			}
		}else{
			diff = 0.;
			printf("Parallel Update: \n");
            
//             v_t=max(abs(A).^2*x_var,1e-12);
//             w_t=A*x_t-g.*v_t;
//             [g,dg] = eta_out(w_t,v_t,y,sigma2);
//             s_t=alpha*s_t+(1-alpha)*max((abs(A').^2*(-dg)).^-1,1e-12);
//             r_t=alpha*r_t+(1-alpha)*(x_t+s_t.*(A'*g));
//             [x_t,x_var]=ComplexBournoulliGaussian( r_t, s_t, prmts );

			//Perform a normal GAMP update here.  Use only the variables a, c, and g
			//Still need to figure out how to remove dependence on dg
			// for (key = 0; key < n; key++) {
				// i = seq[key];
				
			//Equivalent matlab code below
			// v = sqrF * c;
			// w = F * a - v .* g;
			// [g, dg, channel_prmts] = channel(y, w, v, channel_prmts);
				
			/* Update a_proj and c_proj */
			for (mu = 0; mu < m; mu++)
				a_proj[mu] = c_proj[mu] = 0;
			for (i = 0; i < n; i++) {
				for (idx = jc[i]; idx < jc[i + 1]; idx++) {
					a_proj[ ir[idx] ] += F[idx] * a[i];
					c_proj[ ir[idx] ] += (conj(F[idx]) * F[idx]) * c[i];
				}
			}

			/* Update w and v */
			for (mu = 0; mu < m; mu++) {
				o[mu] = a_proj[mu] - c_proj[mu] * g[mu];
				v[mu] = c_proj[mu];
			}
            
            if (t==0){
                printf("At t=0 o[0] is %f + i%f\n", creal(o[0]), cimag(o[0]));
                printf("At t=0 v[0] is %.4e \n", v[0]);
            }
                
            
            for (mu = 0; mu < m; mu++)
				channel(y[mu], o[mu], v[mu], delta, &g[mu], &dg[mu], vnf);

			//for (mu = 0; mu < m; mu++) g_old[mu] = g[mu];
				
			/* Update {sig, r}, {a, c} */
			for (i = 0; i < n; i++){ 
				//a_old_vec[i] = a[i], c_old = c[i];

				g_proj = dg_proj = 0.; /* Dot products: F * g and -F^2 * dg */
				for (idx = jc[i]; idx < jc[i + 1]; idx++) {
					g_proj += F[idx] * g[ ir[idx] ];
					dg_proj -= (conj(F[idx]) * F[idx]) * dg[ ir[idx] ];
				}
                
                dg_proj_vec[i]=dg_proj;
                
// 				if(t>0)
					sig[i] = damp_par * sig[i] + (1 - damp_par) * (1. / dg_proj);
// 				else
// 					sig[i] = 1. / dg_proj;

				if(sig[i] < 0) sig[i] = 0.1*vnf;
				//if(sig[i] < .1*vnf) sig[i] = 0.1*vnf;

// 				if(t>0)
					r[i] = damp_par * r[i] + (1 - damp_par) * (a[i] + sig[i] * g_proj);
// 				else
// 					r[i] = a[i] + sig[i] * g_proj;
			}
            
            if (t==0){
                printf("At t=0 dg_proj_vec[0]^-1 %f + i%f\n", creal(1./dg_proj_vec[0]), cimag(1./dg_proj_vec[0]));
                printf("At t=0 g[0] is %f + i%f\n", creal(g[0]), cimag(g[0]));
                printf("At t=0 dg[0] is %.4e \n", dg[0]);
                printf("At t=0 damp_par is %.4e \n", damp_par);
            }
            
            if (t==0){
                printf("Before update, at t=0 r[0] is %f + i%f\n", creal(g[0]), cimag(g[0]));
                printf("Before update, at t=0 sig[0] is %.4e \n", dg[0]);
                printf("Before update, at t=0 a[0] is %f + i%f\n", creal(a[0]), cimag(a[0]));
                printf("Before update, at t=0 c[0] is %.4e \n", c[0]);
            }
            
			/* Update {a,c} */
/*			for (i = 0; i < n; i++){
				prior_nzmean(r[i], sig[i], p_prms, &a[i], &c[i], vnf);
				//prior(r[i], sig[i], p_prms, &a[i], &c[i], vnf);
				//a[i]=(1.0-.10)*a_old_vec[i]+0.10*a[i];
			}
*/
			prior_nzmean(n, r, sig, p_prms, a, c, vnf);
			
            if (t==0){
                printf("At t=0 r[0] is %f + i%f\n", creal(g[0]), cimag(g[0]));
                printf("At t=0 sig[0] is %.4e \n", dg[0]);
                printf("At t=0 a[0] is %f + i%f\n", creal(a[0]), cimag(a[0]));
                printf("At t=0 c[0] is %.4e \n", c[0]);
            }
			
				// /* Update {w, v}, {g, dg} */
				// for (idx = jc[i]; idx < jc[i + 1]; idx++) {
					// mu = ir[idx];
					// v_old = v[mu];

					// v[mu] += (F[idx] * F[idx]) * (c[i] - c_old);
					// o[mu] += F[idx] * (a[i] - a_old)
							// - (v[mu] - v_old) * g_old[mu];

					// channel(y[mu], o[mu], v[mu], delta, &g[mu], &dg[mu], vnf);
				// }
			for (i = 0; i < n; i++){
				diff += cabs(a[i] - a_old);
			}
		}
        diff /= n;
        
        /* The following helps with debugging. Can compare results to matlab version */
//         channel(3.0, 1.0, 2.0, 4.0, &g_test[1], &dg_test[1], vnf);
//         prior(1.0,2.0, p_prms, &a_test[1], &c_test[1], vnf);
// 
//         printf("a_test[1] is %f + i%f\n", creal(a_test[1]), cimag(a_test[1]));
//         printf("c_test[1] is %.4e \n", c_test[1]);
//         printf("g_test[1] is %.4e \n", g_test[1]);
//         printf("dg_test[1] is %.4e \n", dg_test[1]);
// 		printf("dg[1] is %.4e \n", dg[1]);
		
        dg_ave=0.0;
		for (i = 0; i < n; i++)
			dg_ave+=dg_proj_vec[i];
		dg_ave=dg_ave/n;
		printf("Average dg_proj is %.4e \n", dg_ave);
        
		sig_ave=0.0;
		for (i = 0; i < n; i++)
			sig_ave+=sig[i];
		sig_ave=sig_ave/n;
		printf("Average Sig is %.4e \n", sig_ave);

		r_norm=0.0;
		for (i = 0; i < n; i++)
			r_norm+=pow(cabs(r[i]),2);
		r_norm=sqrt(r_norm);
		printf("Norm of R is %.4e \n", a_norm);
		
		a_norm=0.0;
		for (i = 0; i < n; i++)
			a_norm+=pow(cabs(a[i]),2);
		a_norm=sqrt(a_norm);
		printf("Norm of Estimate is %.4e \n", a_norm);
	
        for (mu = 0; mu < m; mu++)
            a_proj[mu] = 0;
        for (i = 0; i < n; i++) {
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += F[idx] * a[i];
            }
        }
        diff_y = 0.0;
        for (mu = 0; mu < m; mu++){
            diff_y += fabs(y[mu] - cabs(a_proj[mu]));
			// if ((mu%10000)==0){
				// temp=y[mu] - cabs(a_proj[mu]);
				// printf("y[%3d] is %.4e, and abs(a_proj[%3d]) is %.4e, and difference is %.4e \n", mu, y[mu], mu, cabs(a_proj[mu]), temp);
				// printf("Intermediate (iter %3d) diff_y %.4e \n", mu, diff_y);
			// }
		}
        //diff_y /= m;

        if (output)
            fprintf(output, "%g\n", diff_y);
        mexEvalString("drawnow");
		


        /* Print some info. */
        if (disp){
            printf("t: %3d; diff: %.4e; diff_y: %.4e\n \n", t, diff, diff_y);
// 			printf("y[0]: %.4e; |a_proj[0]|: %.4e\n", y[0], cabs(a_proj[0]));
// 			temp=y[0] - cabs(a_proj[0]);
// 			printf("y[0]-abs(a_proj[0]): %.4e\n", temp);
			//printf("Enter ctrl+c (or ctrl+c+j) to Exit\n \n");  
			//myChar=getchar();
			//if (myChar=='\n' || myChar==EOF) break;
		}//mexEvalString("drawnow");
		
        /* Check for convergence */
        //if (diff < eps || isnan(diff)) break;
		

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
	
	free(a_test);
	free(c_test);
	free(g_test);
	free(dg_test);
	free(dg_proj_vec);
	free(a_old_vec);
}
