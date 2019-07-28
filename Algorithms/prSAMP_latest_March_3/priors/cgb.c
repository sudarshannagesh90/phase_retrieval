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


void prior( complex double r, double sig, double *prmts,
        complex double *a, double *c, double vnf ) {
    double rho = prmts[0];
    double xm = prmts[1];
    double xv = prmts[2];
    double z, s, alpha;
    complex double M ;



        s = xv*sig/(sig + xv);
        M = (xv*r + sig * xm) / (sig + xv);

        alpha = pow(cabs(xm),2)/xv - pow(cabs(M),2)/s;

        z = (1-rho)*exp(alpha/2) + rho*sig/(sig+xv);
        if(z<1e-12)  z = 0.1*vnf;
        (*a) = rho*sig/(sig+xv)*M;
        (*a) = (*a)/z;

        (*c) = (1/z)*rho*sig/(sig+xv)*(2*s + pow(cabs(M),2)) - pow(cabs(*a),2);
        (*c) = 0.5* (*c);

        if(*c<1e-18 || !isfinite(*c)) *c = 0.1*vnf;


}


