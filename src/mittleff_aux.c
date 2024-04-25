#include <complex.h>
#include <stdio.h>
#include <math.h>

extern double
compute_r1 (const double alpha, const double acc);

extern double complex
asymptotics (const double alpha, const double beta,
             const double complex z, const double acc);

extern double complex
integral_rep (const double alpha, const double beta,
             const double complex z, const double acc);

double complex 
mittleff_aux (const double alpha,
              const double beta,
              const double complex z,
              const double acc,
              const double r0)
{    
    double complex res;
    double r1;

    r1 = compute_r1(alpha, acc);
    //printf("r1 = %g\n", r1);
    
    if (cabs(z) >= r1) // 1-4        
        res = asymptotics(alpha, beta, z, acc);
    else // 5-6
        res = integral_rep(alpha, beta, z, acc);
    return res;
}
