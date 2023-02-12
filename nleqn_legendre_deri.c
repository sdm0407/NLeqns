#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "C:\Program Files\Simcenter\2021.2\Amesim\lib\ameutils.h"

#define Lgd_Ord 5

static unsigned cnt = 0;

static double Legendre(double x, int n)
{
    double p_0, p_1;
    double temp;

    p_0 = 1.;
    p_1 = x;

    for (int i = 2; i < n + 1; i++)
    {
        temp = ((2 * i - 1) * x * p_1 - (i - 1) * p_0) / i;
        p_0 = p_1;
        p_1 = temp;
    }

    ++cnt;

    return p_1;
}

static double dLegndre(double *x, void *data)
{
    int *n = (int *)data;

    double p_0, p_1;
    double dp_0, dp_1;
    double temp;

    p_0 = 1.;
    p_1 = *x;

    dp_0 = 0.;
    dp_1 = 1;

    for (int i = 2; i < (*n) + 1; i++)
    {
        temp = ((2 * i - 1) * p_1 + (2 * i - 1) * (*x) * dp_1 - (i - 1) * dp_0) / i;
        dp_0 = dp_1;
        dp_1 = temp;

        temp = ((2 * i - 1) * (*x) * p_1 - (i - 1) * p_0) / i;
        p_0 = p_1;
        p_1 = temp;
    }

    ++cnt;

    return dp_1;
}

int main(int argc, char *argv[])
{
    int max_num_it, found;
    double x, xl, xr;
    int userdata[1];

    double roots[Lgd_Ord + 1];
    double Lvals[Lgd_Ord + 1];

    double qq;

    qq = Lgd_Ord * (Lgd_Ord + 1);

    double ws[Lgd_Ord + 1];

    double div1, sub1;

    double mD[Lgd_Ord + 1][Lgd_Ord + 1];

    double vI[Lgd_Ord + 1];

    max_num_it = 200;

    /* polynomial coefficients. */
    userdata[0] = Lgd_Ord;

    unsigned ndiv = 200;
    double dx = 2. / (ndiv + 1.);

    double y0, y1;

    unsigned idx = 1;

    xl = -1 - dx;
    xr = xl + dx;

    y0 = dLegndre(&xr, userdata);

    roots[0] = -1.0;
    roots[Lgd_Ord] = 1.;

    Lvals[0] = Legendre(-1.0, Lgd_Ord);
    Lvals[Lgd_Ord] = Legendre(1.0, Lgd_Ord);

    for (int i = 0; i < ndiv; ++i)
    {
        xl += dx;
        xr = xl + dx;

        y1 = dLegndre(&xr, userdata);

        if (y0 * y1 > 0)
            continue;

        cnt = 0;

        ameNLSolve_(&xl, &x, &xr, dLegndre, &userdata, &max_num_it, &found);

        if (found == 1)
        {
            roots[idx] = x;
            Lvals[idx] = Legendre(x, Lgd_Ord);

            fprintf(stderr, "the %2dth solution = %g.\n iteration = %d\n", idx, x, cnt);
        }
        else
        {
            fprintf(stderr, "Convergence not obtained using %d iterations\n", max_num_it);
        }

        ++idx;

        y0 = y1;
    }

    for (int i = 0; i < Lgd_Ord + 1; ++i)
    {

        ws[i] = 2 / qq / Lvals[i] / Lvals[i];

        for (int j = 0; j < Lgd_Ord + 1; ++j)
        {
            div1 = Lvals[i] / Lvals[j];

            sub1 = roots[i] - roots[j];

            if (i == j)
                sub1 += 1;

            mD[i][j] = div1 / sub1;

            if (i == j)
                mD[i][j] = 0;
        }

        vI[i] = 0.;
    }

    mD[0][0] = -qq / 4;
    mD[Lgd_Ord][Lgd_Ord] = qq / 4;

    vI[0] = vI[Lgd_Ord] = 1.0;

    return 0;
}