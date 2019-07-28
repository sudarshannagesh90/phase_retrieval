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

#include "../arbmat_cswgamp_newchannel.h"


void prior_nzmean( size_t n, complex double *r, double *sig, double *prmts, complex double *a, double *c, double vnf ) {
    double rho = prmts[0];
    //double xm = prmts[1];
    double xv = prmts[2];
	double *xm_array;
	double xm;
    double z, s, alpha;
    complex double M ;
	double rescaling;
	unsigned int i;
	
	xm_array = malloc(sizeof(double) * n);
	
	rescaling=(double)0.0/n;
	for (i = 0; i < n; i++)
		xm_array[i]=(i+1)*rescaling;//should producing float result
	
	for (i = 0; i < n; i++){
		xm=xm_array[i];
        s = xv*sig[i]/(sig[i] + xv);
        M = (xv*r[i] + sig[i] * xm) / (sig[i] + xv);

        alpha = pow(cabs(xm),2)/xv - pow(cabs(M),2)/s;

        z = (1-rho)*exp(alpha/2) + rho*sig[i]/(sig[i]+xv);
        if(z<1e-12)  z = 0.1*vnf;
        (a[i]) = rho*sig[i]/(sig[i]+xv)*M;
        (a[i]) = (a[i])/z;

        (c[i]) = (1/z)*rho*sig[i]/(sig[i]+xv)*(2*s + pow(cabs(M),2)) - pow(cabs(a[i]),2);
        (c[i]) = 0.5* (c[i]);

        if(c[i]<1e-18 || !isfinite(c[i])) c[i] = 0.1*vnf;
	}
	
	free(xm_array);

}


