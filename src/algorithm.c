#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

extern double complex
ml_erfc (double complex x, double eps);

double complex
series (const double alpha,
        const double beta,
        const double complex z,
        const double acc,
        int asymptotic);

extern double
rgamma (const double x, const double eps);
extern int
in_region_G1 (double complex z, double a);
extern int
in_region_G2 (double complex z, double a);
extern int
in_region_G3 (double complex z, double a);
extern int
in_region_G4 (double complex z, double a);
extern int
in_region_G5 (double complex z, double a, double r0);
extern int
in_region_G6 (double complex z, double a, double r0);

double complex
taylor (const double alpha,
        const double beta,
        const double complex z,
        const double acc)
{
    return series(alpha, beta, z, acc, 0);
}

double complex
c_theta (double th)
{
    double complex x, y;

    x = 2.0 * pow(0.5 * th, 2);
    y = th - sin(th);
    
    return csqrt(2.0 * (x + I*y));
}

double complex
rho_i (int i, double complex z, double alpha, double acc)
{
    double complex res;
    double th;

    th = carg(z);
    if (i == 1)
        res = 1.0;
    else if (i == 2)
        res = 0.0;
    else if (i == 3)
        res = 0.5 * ml_erfc(+1 * c_theta(th) * csqrt(0.5 * pow(cabs(z), 1/alpha)), acc);
    else if (i == 4)
        res = 0.5 * ml_erfc(-1 * c_theta(th) * csqrt(0.5 * pow(cabs(z), 1/alpha)), acc);
    return res;
    
}

/* double complex */
/* fn_a (double complex z, */
/*       double a, */
/*       double b, */
/*       double complex x) */
/* { */
/*     /\* double complex tmp = (1.0 / a) * cpow(z, (1 - b)/a) * cexp(cpow(z, 1/a) * ccos(x/a)); *\/ */
/*     /\* printf("x, a = (%g %g)\n", creal(tmp), cimag(tmp));//return 0; *\/ */
/*     return (1.0 / a) * cpow(z, (1 - b)/a) * cexp(cpow(z, 1/a) * ccos(x/a)); */
/* } */

double complex
omega (double complex x,
       double complex y,
       double a,
       double b)
{
    return (cpow(x, 1.0/a) * csin(y/a)) + (y * (1 + ((1.0 - b)/a)));
}

double complex
asymptotics_aux (const double alpha,
                 const double beta,
                 const double complex z,
                 const int region,
                 const double acc)
{
    double complex res, fac1, fac2, fac3;

    fac1 = rho_i(region, z, alpha, acc)/alpha;
    fac2 = cpow(z, (1 - beta)/alpha);
    fac3 = cexp(cpow(z, 1/alpha));
    
    res = fac1*fac2*fac3 - series(alpha, beta, z, acc, 1);

    return res;
}

double complex
asymptotics (const double alpha,
             const double beta,
             const double complex z,
             const double acc)
{
    int region;
    
    if (in_region_G1(z, alpha))
        region = 1;
    else if (in_region_G2(z, alpha))
         region = 2;
    else if (in_region_G3(z, alpha))
         region = 3;
    else if (in_region_G4(z, alpha))
        region = 4;

    assert(region == 1 || region == 2 || region == 3 || region == 4);

    return asymptotics_aux(alpha, beta, z, region, acc);
}

double
Phi_i (int region, double alpha)
{
    double res;
    
    if (region == 5)
        res = M_PI * alpha;
    else if (region == 6)
        res = (2.0/3.0) * M_PI * alpha;
    
    return res;
}

double
Rho_beta (double beta)
{
    return (beta > 1) ? 0.5 : 0.0;
}

double
max3 (double a, double b, double c)
{
    double res;

    res = (a > b)?a:b;
    res = (res > c)?res:c;
    return res;
}

double
Rib (int region, double a, double b,
     double complex z,
    double acc)
{
    double res, r1, r2, r3;

    if ((b >= 0) && (region == 5))
    {
        r1 = 1;
        r2 = 2.0 * cabs(z);
        r3 = pow(-log(M_PI*acc/6.0), a);            
    }
    else if ((b < 0) && (region == 5))
    {
        r1 = pow(1+fabs(b), a);
        r2 = 2.0 * cabs(z);
        r3 = pow(-2.0 * log((M_PI*acc)/(6*(2+fabs(b))*pow(2.0*fabs(b), fabs(b)))), a);
    }
    else if ((b >= 0) && (region == 6))
    {
        r1 = pow(2, a);
        r2 = 2*cabs(z);
        r3 = pow(-2*log(M_PI*acc*(1/12.)*pow(2,b)),a);
    }
    else if ((b < 0) && (region == 6))
    {
        r1 = pow(2*(1+fabs(b)),a);
        r2 = 2*cabs(z);
        r3 = pow(-4*log((M_PI*acc*pow(2,b))/(12*(2+fabs(b))*(pow(4*fabs(b), fabs(b))))),a);
    }
    
    res = max3(r1, r2, r3);

    return res;
}

double complex
integral_rep_aux (const double alpha,
                 const double beta,
                 const double complex z,
                 const int region,
                 const double acc)
{
    return 0;
}

double complex
integral_rep (const double alpha,
              const double beta,
              const double complex z,
              const double acc)
{
    int region;

    if (in_region_G5(z, alpha, 0.95))
        region = 5;
    else if (in_region_G6(z, alpha, 0.95))
        region = 6;
    assert(region == 5 || region == 6);
    
    return integral_rep_aux(alpha, beta, z, region, acc);
}


void
compute_kmin_kmax(int *s, int *kmin, int *kmax,
                  const double alpha, const double beta,
                  const double complex z, const double acc,
                  int asymptotic)
{
    int k1, k2;

    if (asymptotic)
    {
        *s = -1;
        *kmin = 1;        
        *kmax = (int)(ceil((1/alpha) * pow(cabs(z), 1/alpha)) + 1);               
    }
    else
    {
        *s = +1;
        *kmin = 0;

        k1 = (int)(ceil((2+beta)/alpha) + 1);
        k2 = (int)(ceil(log(acc * (1-cabs(z)))/log(cabs(z))) + 1);
        *kmax = (k1 > k2) ? k1 : k2;        
    }
}

double complex
series (const double alpha,
        const double beta,
        const double complex z,
        const double acc,
        int asymptotic)
{
    double complex res, zk;
    double rg;
    int k, kmin, kmax, s;

    compute_kmin_kmax(&s, &kmin, &kmax, alpha, beta, z, acc, asymptotic);

    res = 0;
    for (k = kmin; k <= kmax; k++)
    {
        zk = cpow(z, s*k);
        rg = rgamma(alpha * s * k + beta, acc);

        if (cabs(zk) == 0 || fabs(rg) == 0)
            res += 0;
        else 
            res += zk * rg;
    }
    
    return res;
}


