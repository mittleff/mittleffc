#include <complex.h>
#include <math.h>
#include <stdio.h>

extern double complex 
mittleff_aux (const double alpha,
              const double beta,
              const double complex z,
              const double acc,
              const double r0);

extern double complex
taylor (const double alpha,
        const double beta,
        const double complex z,
        const double acc);

double complex
mittleff (double alpha, double beta, double complex z, double acc, double taylor_radius)
{
    //printf("%g %g %g %g %g %g\n", alpha, beta, creal(z), cimag(z), acc, taylor_radius);
    double r0;
    double complex res;

    r0 = taylor_radius;

    if (cabs(z) <= r0)
    {
        printf("aaa\n");
        // Taylor
        res = taylor(alpha, beta, z, acc);
    }
    else
    {
        if (alpha > 1) // alpha > 1
        {
            // Recursive
            res = 0;
        }
        else
        {
            printf("bbb\n");
            // Main algorithm
            res = mittleff_aux(alpha, beta, z, acc, r0);
        }
    }

    return res;
}
