#include <stdio.h>
#include <math.h>

double P0(double x)
{
    return 1;
}

double P1(double x)
{
    return x;
}
// The following is a general functoin that returns the value of the Legendre Polynomial for any given x and n=0,1,2,3,...
double Pn(double x, int n)
{
    if (n == 0)
    {
        return P0(x);
    }
    else if (n == 1)
    {
        return P1(x);
    }
    else
    {
        return (double)((2 * n - 1) * x * Pn(x, n - 1) - (n - 1) * Pn(x, n - 2)) / n;
    }
}

double myPn(double x, int n)
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

    return p_1;
}

double mydPn(double x, int n)
{
    double p_0, p_1;
    double dp_0, dp_1;
    double temp;

    p_0 = 1.;
    p_1 = x;

    dp_0 = 0.;
    dp_1 = 1;

    for (int i = 2; i < n + 1; i++)
    {
        temp = ((2 * i - 1) * p_1 + (2 * i - 1) * x * dp_1 - (i - 1) * dp_0) / i;
        dp_0 = dp_1;
        dp_1 = temp;

        temp = ((2 * i - 1) * x * p_1 - (i - 1) * p_0) / i;
        p_0 = p_1;
        p_1 = temp;
    }

    return dp_1;
}

int main()
{
    // We will create a data-file and store the values of first few Legendre polynomials for -1<x<1
    FILE *fp1 = NULL, *fp2 = NULL, *fp3 = NULL;
    // create data-file
    fp1 = fopen("legendre1.txt", "w");
    double x;
    // write the values of first 5 Legendre polynomials to data-file
    for (x = -1; x <= 1; x = x + 0.1)
    {
        fprintf(fp1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, Pn(x, 1), Pn(x, 2), Pn(x, 3), Pn(x, 4), Pn(x, 5));
    }
    fclose(fp1);

    fp2 = fopen("legendre2.txt", "w");

    // write the values of first 5 Legendre polynomials to data-file
    for (x = -1; x <= 1; x = x + 0.1)
    {
        fprintf(fp2, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, myPn(x, 1), myPn(x, 2), myPn(x, 3), myPn(x, 4), myPn(x, 5));
    }
    fclose(fp2);

    fp3 = fopen("legendre2_deri.txt", "w");

    // write the values of first 5 Legendre polynomials to data-file
    for (x = -1; x <= 1; x = x + 0.1)
    {
        fprintf(fp3, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, mydPn(x, 1), mydPn(x, 2), mydPn(x, 3), mydPn(x, 4), mydPn(x, 5));
    }

    fclose(fp3);

    return 0;
}