#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "C:\Program Files\Simcenter\2021.2\Amesim\lib\ameutils.h"

typedef struct
{
    int sz;
    double *roots;
    double *mD;
    double *vI;
    double *ws;
    double *mD2;
} SEMdata;

static unsigned cnt = 0;

double Legendre(double, int);

double dLdr(double, int);

SEMdata *SEM2(int);

int bisect(double, double, double (*dLdr)(double, int), int, double *);

static void fun15(double *x, int *n, void *data, double *u, int *ok);

static unsigned f_iter, j_iter;

#define ndiv 200 // number of division for bisecting the domain
#define tol 1e-7

#define Origin -1
#define nZones 2

#define Lgd_Ord 4
#define nElem 4

const int nO[nZones] = {4, 5};
const int nE[nZones] = {2, 2};

const double Leng[nZones] = {1.0, 1.0};

int main(int argc, char *argv[])
{

    int jNO1[nZones];

    for (int i = 0; i < nZones; ++i)
    {
        jNO1[i] = nO[i] * nE[i];
    }

    double dlV[nZones];
    double MxdlV = 0.0;
    for (int i = 0; i < nZones; ++i)
    {
        dlV[i] = Leng[i] / nE[i];

        if (dlV[i] > MxdlV)
            MxdlV = dlV[i];
    }

    int nU = 0;
    double TotLeng;
    for (int i = 0; i < nZones; ++i)
    {
        nU += jNO1[i];
        TotLeng += Leng[i];
    }

    nU = nU + 1;

    SEMdata *mSE[nZones];
    for (int i = 0; i < nZones; ++i)
    {
        mSE[i] = SEM2(nO[i]);
    }

    // # Calculate the solutions: Y value

    double tolf, tolx;
    int max_num_it, found, num_steps;

    tolf = 1.0e-6;
    tolx = 1.0e-7;
    max_num_it = 100;

    double *yy = (double *)malloc(sizeof(double) * nU);
    if (yy == NULL)
    {
        fprintf(stderr, "x can't be allocated\n");
        return 1;
    }

    for (int i; i < nU; ++i)
        yy[i] = 1.0;

    num_steps = ndnr_(yy, &nU, fun15, NULL, mSE, &tolf, &tolx, &max_num_it, &found, NULL);

    // # Calculate the domain position: X value
    double *xx = (double *)malloc(sizeof(double) * nU);
    if (xx == NULL)
    {
        fprintf(stderr, "x can't be allocated\n");
        return 1;
    }

    int idx0 = 0, idx1 = 0;

    double x0 = Origin;
    double e_l;

    int sz;

    for (int i = 0; i < nZones; ++i)
    {
        sz = mSE[i]->sz;
        double *dx = (double *)alloca(sizeof(double) * sz);

        e_l = (Leng[i] / nE[i]) / 2;

        for (int j = 0; j < sz; ++j)
        {
            dx[j] = e_l * (mSE[i]->roots[j] - 1) + x0;
        }

        for (int k = 0; k < nE[i]; ++k)
        {
            idx0 = idx1;
            idx1 += nO[i];

            for (int j = 0; j < sz; ++j)
            {
                dx[j] += 2. * e_l;
                xx[idx0 + j] = dx[j];
            }
        }

        x0 = dx[sz - 1];
    }

    // # Calculate the solutions: dY value
    double *dy = (double *)malloc(sizeof(double) * nU);
    if (dy == NULL)
    {
        fprintf(stderr, "dy can't be allocated\n");
        return 1;
    }

    double temp;
    idx0 = 0;
    idx1 = 0;

    for (int i = 0; i < nZones; ++i)
    {
        sz = mSE[i]->sz;

        for (int k = 0; k < nE[i]; ++k)
        {
            idx0 = idx1;
            idx1 += nO[i];

            for (int i1 = 0; i1 < sz; ++i1)
            {
                temp = 0.;
                for (int j1 = 0; j1 < sz; ++j1)
                {
                    temp += mSE[i]->mD[i1 * sz + j1] * yy[idx0 + j1];
                }
                dy[idx0 + i1] = temp / e_l;
            }
        }
    }

    // print results
    FILE *xfptr = fopen("results\\x.csv", "w");
    FILE *yfptr = fopen("results\\y.csv", "w");
    FILE *dyfptr = fopen("results\\dy.csv", "w");

    for (int i = 0; i < nU; ++i)
    {
        fprintf(xfptr, "%f\n", xx[i]);
        fprintf(yfptr, "%f\n", yy[i]);
        fprintf(dyfptr, "%f\n", dy[i]);
    }

    fclose(xfptr);
    fclose(yfptr);
    fclose(dyfptr);

    printf("iteration(f_iter) is %d\n", f_iter);
    // free all the memory allocated memory in Heap before exit the program.
    for (int i = 0; i < nZones; ++i)
    {
        free(mSE[i]->roots);
        free(mSE[i]->mD);
        free(mSE[i]->vI);
        free(mSE[i]->ws);
        free(mSE[i]->mD2);
        free(mSE[i]);
        mSE[i] = NULL;
    }

    free(xx);
    xx = NULL;

    free(yy);
    yy = NULL;

    free(dy);
    dy = NULL;

    return 0;
}

double Legendre(double x, int n)
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

double dLdr(double x, int n)
{

    double p0, p1;
    double dp0, dp1;
    double temp;

    p0 = 1.;
    p1 = x;

    dp0 = 0.;
    dp1 = 1;

    for (int i = 2; i < n + 1; i++)
    {
        temp = ((2 * i - 1) * p1 + (2 * i - 1) * x * dp1 - (i - 1) * dp0) / i;
        dp0 = dp1;
        dp1 = temp;

        temp = ((2 * i - 1) * x * p1 - (i - 1) * p0) / i;
        p0 = p1;
        p1 = temp;
    }

    ++cnt;

    return dp1;
}

SEMdata *SEM2(int Ord)
{

    double x, xl, xr;

    int ord1 = Ord + 1;

    SEMdata *p = (SEMdata *)malloc(sizeof(SEMdata));
    if (p == NULL)
        return NULL;

    p->sz = ord1;

    p->roots = (double *)malloc(sizeof(double) * ord1);
    if (p->roots == NULL)
        return NULL;

    p->mD = (double *)malloc(sizeof(double) * ord1 * ord1);
    if (p->mD == NULL)
        return NULL;

    p->vI = (double *)malloc(sizeof(double) * ord1);
    if (p->vI == NULL)
        return NULL;

    p->ws = (double *)malloc(sizeof(double) * ord1);
    if (p->ws == NULL)
        return NULL;

    p->mD2 = (double *)malloc(sizeof(double) * ord1 * ord1);
    if (p->mD2 == NULL)
        return NULL;

    double *Lvals = (double *)alloca(sizeof(double) * ord1); // on Stack
    if (Lvals == NULL)
        return NULL;

    double qq = Ord * ord1;

    double div1, sub1;

    /* polynomial coefficients. */

    double dx = 2. / (ndiv - 1.);

    double y0, y1;

    unsigned idx = 1;

    xl = -1 - dx;
    xr = xl + dx;

    y0 = dLdr(xr, Ord);

    p->roots[0] = -1.0;
    p->roots[Ord] = 1.0;

    Lvals[0] = Legendre(-1.0, Ord);
    Lvals[Ord] = Legendre(1.0, Ord);

    for (int i = 0; i < ndiv; ++i)
    {
        xl += dx;
        xr = xl + dx;

        y1 = dLdr(xr, Ord);

        if (y0 * y1 > 0)
            continue;

        cnt = 0;

        cnt = bisect(xl, xr, dLdr, Ord, &x);

        if (cnt)
        {
            p->roots[idx] = x;
            Lvals[idx] = Legendre(x, Ord);

            printf("the %2dth solution = %g. iteration = %d\n", idx, x, cnt);
        }
        else
        {
            printf("Convergence not obtained \n");
        }

        ++idx;

        y0 = y1;
    }

    for (int i = 0; i < ord1; ++i)
    {

        p->ws[i] = 2 / qq / Lvals[i] / Lvals[i];

        for (int j = 0; j < ord1; ++j)
        {
            div1 = Lvals[i] / Lvals[j];

            sub1 = p->roots[i] - p->roots[j];

            if (i == j)
                sub1 += 1;

            p->mD[i * ord1 + j] = div1 / sub1;

            if (i == j)
                p->mD[i * ord1 + j] = 0;
        }

        p->vI[i] = 0.;
    }

    p->mD[0] = -qq / 4;
    p->mD[ord1 * ord1 - 1] = qq / 4;

    p->vI[0] = p->vI[Ord] = 1.0;

    for (int i = 0; i < ord1; ++i)
    {

        for (int j = 0; j < ord1; ++j)
        {
            p->mD2[i * ord1 + j] = 0.0;

            for (int k = 0; k < ord1; ++k)
                p->mD2[i * ord1 + j] += p->mD[i * ord1 + k] * p->mD[j + k * ord1];
        }
    }

    return p;
}

int bisect(double xmin, double xmax, double (*dLdr)(double, int), int n, double *sol)
{
    double x;

    double y0, y1;

    y0 = dLdr(xmin, n);

    int cnt1 = 0;

    while (true)
    {
        x = (xmin + xmax) / 2.;

        y1 = dLdr(x, n);

        if (y0 * y1 > 0)
        {
            xmin = x;
            y0 = y1;
        }

        else
            xmax = x;

        ++cnt1;

        if (fabs(y1) < tol)
            break;
    }

    *sol = x;

    return cnt1;
}

/* Computes F(x) */

static void fun15(double *x, int *n, void *data, double *u, int *ok)
{
    SEMdata **sem;

    sem = (SEMdata **)data;

    int idx0 = 0, idx1 = 0;

    double x0 = Origin;

    int sz;

    for (int i = 0; i < *n; ++i)
        u[i] = 0.0;

    // checking values
    // double m1D[Lgd_Ord + 1][Lgd_Ord + 1];
    // double m2D[Lgd_Ord + 1][Lgd_Ord + 1];
    // double ys[Lgd_Ord + 1];

    for (int i = 0; i < nZones; ++i) // loop for zones
    {
        sz = sem[i]->sz;

        double *dx = (double *)alloca(sizeof(double) * sz);

        double e_l = (Leng[i] / nE[i]) / 2;

        for (int j = 0; j < sz; ++j)
        {
            dx[j] = e_l * (sem[i]->roots[j] - 1) + x0;
        }

        for (int k = 0; k < nE[i]; ++k) // loop for elements
        {
            idx0 = idx1;
            idx1 += nO[i];

            for (int j = 0; j < sz; ++j)
            {
                dx[j] += 2. * e_l;
            }

            double temp1 = 0.;

            int i1 = 0;

            temp1 = 0.;
            for (int j1 = 0; j1 < sz; ++j1)
            {
                temp1 += sem[i]->mD[i1 * sz + j1] * x[idx0 + j1];
            }
            u[idx0 + i1] += temp1 / e_l;

            for (i1 = 1; i1 < sz - 1; ++i1)
            {
                temp1 = 0.;
                for (int j1 = 0; j1 < sz; ++j1)
                {
                    temp1 += sem[i]->mD2[i1 * sz + j1] * x[idx0 + j1];
                }
                u[idx0 + i1] = temp1 - exp(4. * dx[i1]) * e_l * e_l;
            }

            temp1 = 0.;
            for (int j1 = 0; j1 < sz; ++j1)
            {
                temp1 += sem[i]->mD[i1 * sz + j1] * x[idx0 + j1];
            }
            u[idx0 + i1] -= temp1 / e_l;
        }

        x0 = dx[sz - 1];
    }

    u[*n - 1] *= -1.0;

    //(* Boundary Conditions *)
    u[0] = x[0] - 0.0; // volt;
    // u[n-1] -= Kls[n-1] * 10   # OutVoltage  # volt
    u[*n - 1] = x[*n - 1] - 0.0; // volt

    *ok = 1;

    ++f_iter;
}
