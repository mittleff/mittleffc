#include <math.h>

#include "algorithm.h"
#include "integrate.h"
#include "log.h"
#include "new.h"

#include <gsl/gsl_math.h>
#include <stdbool.h>

#define GSL_EPSILON_FCMP 1e-10
#define MAX2(a,b) ((a>b)?(a):(b))
#define MAX3(a,b,c) MAX2(a,MAX2(b,c))

static void
asymptotic_series (num_t res, const double x, const double y, const double alpha, const double beta);


static bool
ge (const double x, const double y)
{
    int ret = gsl_fcmp(x, y, GSL_EPSILON_FCMP);
    return (ret == 1 || ret == 0) ? true : false;
}

static bool
le (const double x, const double y)
{
    int ret = gsl_fcmp(x, y, GSL_EPSILON_FCMP);
    return (ret == -1 || ret == 0) ? true : false;
}

void
mittleff0 (double* res,
           const double a, const double b,
           const double x, const double y,
           const double acc)
{
    int k, k1, k2, kmax;
    double abs_z;
    num_t z, sum, tmp, fac1, fac2;

    /* Compute the maximum number of terms kmax to be taken into account for the
     * Taylor series, eq. (4.5) */
    abs_z = sqrt(x*x + y*y);
    k1 = (int) (ceil((2 - b)/a) + 1);
    k2 = (int) (ceil(log(acc * (1 - abs_z))/log(abs_z)) + 1);
    kmax = (k1 > k2) ? k1 : k2;

    /* Sum Taylor series */
    z = new(num), sum = new(num);
    tmp = new(num), fac1 = new(num), fac2 = new(num);
    num_set_d_d(z, x, y), num_set_d(sum, 0.0);
    for (k = 0; k <= kmax; k++)
    {
        /* fac1 <- z**k */
        num_pow_d(fac1, z, (double) k); 

        /* fac2 <- rgamma(alpha * k + beta) */
        num_set_d(fac2, a * k + b); 
        num_rgamma(fac2, fac2);

        /* partial sum = fac1 * fac2 */
        num_mul(tmp, fac1, fac2);
        
        num_add(sum, sum, tmp);
    }
    num_to_d_d(res, sum);
    delete(sum), delete(z), delete(tmp), delete(fac1), delete(fac2);
}

void
mittleff1 (double* res,
           const double alpha, const double beta,
           const double x, const double y,
           const double acc)
{
    /* compute eq. (2.4) */
    num_t z, fac1, fac2, _res;

    z = new(num), fac1 = new(num), fac2 = new(num), _res = new(num);
    num_set_d_d(z, x, y);
    num_pow_d(fac1, z, (1.0 - beta)/alpha);
    num_pow_d(fac2, z, 1.0/alpha);
    num_exp(fac2, fac2);
    num_mul(fac1, fac1, fac2);
    num_mul_d(fac1, fac1, 1.0/alpha);
    asymptotic_series(fac2, x, y, alpha, beta);
    num_add(_res, fac1, fac2);
    num_to_d_d(res, _res);
    delete(_res), delete(fac1), delete(fac2), delete(z);
}

/* num_t */
/* mittleff2 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*     return asymptotic_series(z, alpha, beta); */
/* } */

/* num_t */
/* mittleff3 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*     /\* num_t a = PARAMS_GET_ALPHA(p); *\/ */
/*     /\* num_t b = PARAMS_GET_BETA(p); *\/ */

/*     /\* th = mp.arg(z**(1/alpha)) - np.pi *\/ */
/*     const double th = num_to_d( */
/*         num_sub( */
/*             num_arg(num_pow(z, num_inv(alpha))), */
/*             new(num, M_PI, 0.0))); */
		
/*     /\* c = th + 1j*th**2/6 - th**3/36 *\/ */
/*     const num_t c = new(num, th - th*th*th/36, th*th/6, 0.0); */
		
/*     /\* fac = (1/(2*alpha))*z**((1 - beta)/alpha)*mp.exp(z**(1/alpha))*mp.erfc(c*mp.sqrt(0.5*mp.fabs(z)**(1/alpha))) *\/ */
/*     num_t fac; */
/*     fac = num_inv(num_mul(new(num, 2.0, 0.0), alpha)); */
/*     fac = num_mul(fac, */
/*                   num_pow( */
/*                       z, */
/*                       num_div(num_sub(new(num, 1.0, 0.0), beta), alpha))); */
/*     fac = num_mul(fac, */
/*                   num_exp( */
/*                       num_pow( */
/*                           z, */
/*                           num_inv(alpha)))); */
/*     fac = num_mul(fac, */
/*                   num_erfc( */
/*                       num_mul( */
/*                           c, */
/*                           num_sqrt( */
/*                               num_mul( */
/*                                   new(num, 0.5, 0.0), */
/*                                   num_pow( */
/*                                       num_abs(z), */
/*                                       num_inv(alpha))))))); */
		
/*     return num_add(fac, asymptotic_series(z, alpha, beta)); */
/* } */

/* num_t */
/* mittleff4 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*       /\* numeric_t a = PARAMS_GET_ALPHA(p); *\/ */
/*     /\* numeric_t b = PARAMS_GET_BETA(p); *\/ */

/*     /\* th = mp.arg(z**(1/alpha)) - np.pi *\/ */
/*     const double th = num_to_d( */
/*         num_sub( */
/*             num_arg( */
/*                 num_pow( */
/*                     z, */
/*                     num_inv(alpha))), */
/*             new(num, M_PI, 0.0))); */
		
/*     /\* c = th + 1j*th**2/6 - th**3/36 *\/ */
/*     const num_t c = new(num, -th + th*th*th/36.0, -th*th/6.0); */
		
/*     /\* fac = (1/(2*alpha))*z**((1 - beta)/alpha)*mp.exp(z**(1/alpha))*mp.erfc(c*mp.sqrt(0.5*mp.fabs(z)**(1/alpha))) *\/ */
/*     num_t fac; */
/*     fac = num_inv(num_mul(new(num, 2.0, 0.0), alpha)); */
/*     fac = num_mul( */
/*         fac, */
/*         num_pow( */
/*             z, */
/*             num_div( */
/*                 num_sub(new(num, 1.0, 0.0), beta), */
/*                 alpha))); */
/*     fac = num_mul( */
/*         fac, */
/*         num_exp( */
/*             num_pow( */
/*                 z, */
/*                 num_inv(alpha)))); */
/*     fac = num_mul( */
/*         fac, */
/*         num_erfc(num_mul( */
/*                      c, */
/*                      num_sqrt( */
/*                          num_mul( */
/*                              new(num, 0.5, 0.0), */
/*                              num_pow( */
/*                                  num_abs(z), num_inv(alpha))))))); */
		
/*     return num_add(fac, asymptotic_series(z, alpha, beta));	 */
/* } */

void
mittleff5 (double* res,
           const double alpha, const double beta,
           const double x, const double y,
           const double acc)
{
    log_trace("[%s] (alpha, beta, z, acc) = %+.5e %+.5e (%+.5e%+.5ej) %.5e", __func__, alpha, beta, x, y, acc);
    const double eps = 0.5;

    /* Compute r_max, equation (4.53) */
    double rmax = 0.0;
    double abs_z = sqrt(x*x + y*y);
    if (ge(beta, 0.0))
        rmax = MAX3(1.0,
                    pow(2.0, abs_z),
                    pow(-log(M_PI * eps/6.0), alpha));
    else
        rmax = MAX3(pow(fabs(beta) + 1.0, alpha),
                    2.0*abs_z,
                    pow(-2.0 * log(M_PI*eps/(6.0*(fabs(beta)+2.0)*pow(2.0*fabs(beta), beta))), alpha));

    num_t _res = new(num), a = new(num);    
    A(a, x, y, alpha, beta, 0.0);
    if (le(beta, 1.0)) /* Equation (4.25) */
    {
        num_t integ_b = new(num);
        integrate_B(integ_b, alpha, beta, x, y, M_PI * alpha, 0.0, rmax);
        num_add(_res, a, integ_b);
        delete(integ_b);
    }
    else /* Equation (4.26) */
    {
        num_t integ_b = new(num), integ_c = new(num);
        integrate_B(integ_b, alpha, beta, x, y, M_PI * alpha, 0.5, rmax);
        integrate_C(integ_c, alpha, beta, x, y, 0.5, -M_PI * alpha, M_PI * alpha);
        num_add(_res, integ_b, integ_c);
        num_add(_res, _res, a);
        delete(integ_b), delete(integ_c);
    }
    num_to_d_d(res, _res);
    delete(_res), delete(a);
}

void
mittleff6 (double* res,
           const double alpha,
           const double beta,
           const double x, const double y,
           const double acc)
{
    log_trace("[%s] (alpha, beta, z, acc) = %+.5e %+.5e (%+.5e%+.5ej) %.5e", __func__, alpha, beta, x, y, acc);
    const double eps = 0.5;
    /* Compute r_max, equation (4.54) */
    double rmax = 0.0;
    double abs_z = sqrt(x*x + y*y);
    if (ge(beta, 0.0))
    {
        rmax = MAX3(pow(2.0, alpha),
                    2.0 * abs_z,
                    pow(-log(M_PI*eps*pow(2.0, beta)/12.0), alpha));
    }
    else // TODO Check with Hansjoerg whether the brackets mean "ceil" here
        rmax = MAX3(
            pow(ceil(2.0 * (fabs(beta) + 1.0)), alpha), // HERE
            2.0 * abs_z,
            pow(ceil(-log((M_PI * pow(2.0, beta)*eps)/(12*(fabs(beta) + 2)*pow(4 * fabs(beta), fabs(beta))))), alpha));

    num_t _res = new(num);
    if (le(beta, 1.0)) /* Equation (4.31) */
    {
        num_t integ_b = new(num);
        integrate_B(integ_b, alpha, beta, x, y, 2.0 * M_PI * alpha/3.0, 0.0, rmax);
        num_set(_res, integ_b);
        delete(integ_b);
    }
    else /* Equation (4.32) */
    {
        num_t integ_b = new(num), integ_c = new(num);
        integrate_B(integ_b, alpha, beta, x, y, 2.0 * M_PI * alpha/3.0, 0.5, rmax);
        integrate_C(integ_c, alpha, beta, x, y, 0.5, -2.0 * M_PI * alpha/3.0, 2.0 * M_PI * alpha/3.0);
        num_add(_res, integ_b, integ_c);
        delete(integ_b), delete(integ_c);
    }
    num_to_d_d(res, _res);
    delete(_res);
}


/* rhs of equation (2.3) */
void
asymptotic_series (num_t res, const double x, const double y, const double alpha, const double beta)
{
    int k, kmax;
    double abs_z;
    num_t z, sum, tmp, fac1, fac2;

    abs_z = sqrt(x*x + y*y);
    kmax = (int) (ceil((1.0/alpha) * pow(abs_z, 1.0/alpha)) + 1.0);

    /* -sum([z**(-k) * mp.rgamma(beta - alpha*k) for k in range(1, kmax + 1)]) */
    z = new(num), sum = new(num);
    tmp = new(num), fac1 = new(num), fac2 = new(num);
    num_set_d_d(z, x, y), num_set_d(sum, 0.0);
    for (k = 1; k <= kmax; k++)
    {
        /* fac1 <- z**(-k) */
        num_pow_d(fac1, z, (double) (-k)); 

        /* fac2 <- rgamma(beta - alpha * k) */
        num_set_d(fac2, beta - alpha * k); 
        num_rgamma(fac2, fac2);

        /* partial sum = fac1 * fac2 */
        num_mul(tmp, fac1, fac2);
        
        num_add(sum, sum, tmp);
    }
    num_set(res, sum);
    num_mul_d(res, res, -1.0);
    delete(sum), delete(z), delete(tmp), delete(fac1), delete(fac2);
}
