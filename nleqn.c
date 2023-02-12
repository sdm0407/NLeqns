#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "C:\Program Files\Simcenter\2021.2\Amesim\lib\ameutils.h"

static double userfct(double *x, void *data)
{
    double *coeff = (double *)data;

    return sqrt(*x) - coeff[0] * (*x) * (*x) + coeff[1];
}

int main(int argc, char *argv[])
{
    int max_num_it, found;
    double x, xl, xr, userdata[2];

    max_num_it = 200;

    /* polynomial coefficients. */
    userdata[0] = 3.0;
    userdata[1] = 4.0;

    xl = 0.0;
    xr = 1.e4;

    ameNLSolve_(&xl, &x, &xr, userfct, &userdata, &max_num_it, &found);

    if (found == 1)
    {
        fprintf(stderr, "Result = %g.\n", x);
    }
    else
    {
        fprintf(stderr, "Convergence not obtained using %d iterations\n", max_num_it);
    }

    return 0;
}