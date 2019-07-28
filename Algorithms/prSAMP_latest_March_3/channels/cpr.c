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

void channel(double y, complex double o, double v, double delta,
        complex double *g, double *dg, double vnf) {
    double phi, R0;
    double var;

    if(cabs(o)<1e-12) o = 1e-12;

    phi = 2 * y * cabs(o) / (delta + v);
    R0 = gsl_sf_bessel_I1_scaled(phi) / gsl_sf_bessel_I0_scaled(phi);
    if(phi == 0) R0 = 1;

    *g = o/(delta+v)*(y/cabs(o)*R0-1);

    var = pow(y,2)*(1-R0*R0)/((1+delta/ v)*(1+delta/ v)) + delta * v /(delta+v);
    if(var<1e-12) var = vnf*0.1;
    *dg = 1/ v *(var/ v -1);

}
