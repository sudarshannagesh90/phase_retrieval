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

/* Random order */
void sort_rand( int n, int *seq ) {
    int i, key, exchg;

    for (key = 0; key < n - 1; key++) {
        exchg = key + rand_int(n - key);
        i = seq[key]; seq[key] = seq[exchg]; seq[exchg] = i;
    }
}

double drand(double dMin, double dMax)
{
    return dMin + ( dMax - dMin) * (rand() % RAND_MAX) / RAND_MAX;
}

double logsig(double x)
{
    return 1/(1+exp(-x));
}

//b=A*x  or b=A'*x
//depending on m, n (size of A) and p (size of x)
void project(double *A, double *x, double *b, int m, int n, int p)
{
    int i,j,ind;

   // for(i=0;i<m*n;i++)
   //     printf("index: %d, %.4e\n", i, A[i]);

    if(p==n)
    {
        for (i = 0; i < m; i++) b[i] = 0;
        for (i = 0; i < m; i++)//row
            for (j = 0; j < n; j++) {//column
                ind = i + j * m;
                //printf("index: %d, %.4e\n", ind, A[ind]);
                b[i] += A[ind] * x[j];
            }
    }
    else
    {
        for (i = 0; i < n; i++) b[i] = 0;
        for (i = 0; i < n; i++)//column
            for (j = 0; j < m; j++) {//row
                ind = j + i * m;
                b[i] += A[ind] * x[j];
            }
    }

}
