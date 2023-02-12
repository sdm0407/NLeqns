#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "C:\Program Files\Simcenter\2021.2\Amesim\lib\ameutils.h"

static unsigned f_iter, j_iter;
/* Computes F(x) */

static void userfct(double *x, int *n, void *data, double *u, int *ok)
{
    double *coeff;

    coeff = (double *)data;

    u[0] = coeff[0] - coeff[1] * x[0];
    u[1] = coeff[2] + coeff[3] * (x[1] - coeff[4] * x[0] * x[0]);

    *ok = 1;

    ++f_iter;
}

/* Computes dF(x)/dx jacobian matrix. */

static void userjac(double *x, int *n, void *data, double *jvector, int *ok)
{
    double *coeff;

    coeff = (double *)data;

    jvector[0] = -coeff[1];                         /* df0 / dx0 */
    jvector[1] = -2.0 * coeff[3] * coeff[4] * x[0]; /* df1 / dx0 */
    jvector[2] = 0.0;                               /* df0 / dx1 */
    jvector[3] = coeff[3];                          /* df1 / dx1 */

    *ok = 1;

    ++j_iter;
}

static void userfct2(double *x, int *n, void *data, double *u, int *ok)
{
    double *coeff;

    coeff = (double *)data;
    u[0] = 2. * x[0] + x[1] - 1;
    u[1] = x[0] * x[0] + x[1] * x[1] - 1;

    *ok = 1;

    ++f_iter;
}

static void userjac2(double *x, int *n, void *data, double *jvector, int *ok)
{
    double *coeff;

    coeff = (double *)data;

    jvector[0] = 2.;        /* df0 / dx0 */
    jvector[1] = 2. * x[0]; /* df1 / dx0 */
    jvector[2] = 1.;        /* df0 / dx1 */
    jvector[3] = 2. * x[1]; /* df1 / dx1 */

    *ok = 1;

    ++j_iter;
}

int main(int argc, char *argv[])
{
    int n, max_num_it, found, num_steps;
    double tolf, tolx, x[2], userdata[5];

    tolf = 1.0e-6;
    tolx = 1.0e-7;
    max_num_it = 100;
    n = 2;

    /* polynomial coefficients. */

    userdata[0] = 0.0;
    userdata[1] = 0.7;
    userdata[2] = 1.2;
    userdata[3] = 0.2;
    userdata[4] = -1.0;

    /* Starting point. */

    x[0] = 1.0;
    x[1] = 1.0;

    printf("Start point = (%g, %g)\n", x[0], x[1]);

    x[0] = -1.0;
    x[1] = -1.0;

    f_iter = j_iter = 0;

    num_steps = ndnr_(x, &n, userfct, NULL, &userdata, &tolf, &tolx, &max_num_it, &found, NULL);

    printf("Result = (%g, %g) in %d iterations.\n", x[0], x[1], num_steps);
    printf("Result = f_iter = %d, j_iter = %d \n", f_iter, j_iter);

    x[0] = -1.0;
    x[1] = -1.0;

    f_iter = j_iter = 0;

    num_steps = ndnr_(x, &n, userfct, userjac, &userdata, &tolf, &tolx, &max_num_it, &found, NULL);

    printf("Result = (%g, %g) in %d iterations.\n", x[0], x[1], num_steps);
    printf("Result = f_iter = %d, j_iter = %d \n", f_iter, j_iter);

    x[0] = -1.0;
    x[1] = -1.0;

    f_iter = j_iter = 0;

    num_steps = ndnr_(x, &n, userfct2, NULL, &userdata, &tolf, &tolx, &max_num_it, &found, NULL);

    printf("Result = (%g, %g) in %d iterations.\n", x[0], x[1], num_steps);
    printf("Result = f_iter = %d, j_iter = %d \n", f_iter, j_iter);

    x[0] = -1.0;
    x[1] = -1.0;

    f_iter = j_iter = 0;

    num_steps = ndnr_(x, &n, userfct2, userjac2, &userdata, &tolf, &tolx, &max_num_it, &found, NULL);

    printf("Result = (%g, %g) in %d iterations.\n", x[0], x[1], num_steps);
    printf("Result = f_iter = %d, j_iter = %d \n", f_iter, j_iter);

    system("pause");

    return 0;
}