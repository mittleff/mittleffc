#include <math.h>

#include "algorithm.h"
#include "integrate.h"
#include "log.h"
#include "new.h"

#include <stdbool.h>

#define MAX2(a,b) ((a>b)?(a):(b))
#define MAX3(a,b,c) MAX2(a,MAX2(b,c))

static void
asymptotic_series (num_t res, const num_t z, const num_t alpha, const num_t beta);

void
mittleff0 (num_t res,
           const num_t alpha, const num_t beta,
           const num_t z,
           const num_t acc)
{
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    
    int k, kmax;
    num_t absz, sum, tmp, fac1, fac2, k1, k2;

    absz = new(num);
    k1 = new(num), k2 = new(num);
    fac1 = new(num), fac2 = new(num);
    tmp = new(num);
    sum = new(num);

    /* Compute the maximum number of terms kmax to be taken into account for the
     * Taylor series, eq. (4.5) */
    num_abs(absz, z);
    
    num_one(k1);
    num_set_d(tmp, 2.0);
    num_sub(tmp, tmp, beta);
    num_div(tmp, tmp, alpha);
    num_ceil(tmp, tmp);
    num_add(k1, tmp, k1);
    
    num_one(k2);
    num_one(fac1);
    num_sub(fac1, fac1, absz);
    num_mul(fac1, fac1, acc);
    num_log(fac1, fac1);
    num_log(fac2, absz);
    num_div(tmp, fac1, fac2);
    num_ceil(tmp, tmp);
    num_add(k2, tmp, k2);
    
    //k1 = (int) (ceil((2 - b)/a) + 1);
    //k2 = (int) (ceil(log(acc * (1 - abs_z))/log(abs_z)) + 1);
    kmax = num_gt(k1, k2) ? (int) num_to_d(k1) : (int) num_to_d(k2);
    delete(k1), delete(k2);

    /* Sum Taylor series */
    for (k = 0; k <= kmax; k++)
    {
        /* fac1 <- z**k */
        num_pow_d(fac1, z, (double) k); 

        /* fac2 <- rgamma(alpha * k + beta) */
        num_set_d(fac2, (double) k);
        num_mul(fac2, fac2, alpha);
        num_add(fac2, fac2, beta);
        num_rgamma(fac2, fac2);

        /* partial sum = fac1 * fac2 */
        num_mul(tmp, fac1, fac2);
        
        num_add(sum, sum, tmp);
    }
    num_set(res, sum);
    delete(sum), delete(z), delete(tmp), delete(fac1), delete(fac2);

    //num_print(res, true);
    //log_trace("[%s] Done.", __func__);
    
}

/* compute eq. (2.4) */
void
mittleff1 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__,
          num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc));
    
    num_t fac1, fac2, _res;

    fac1 = new(num), fac2 = new(num), _res = new(num);
    num_one(fac1);
    num_sub(fac1, fac1, beta);
    num_div(fac1, fac1, alpha);
    num_pow(fac1, z, fac1);
    num_inv(fac2, alpha);
    num_pow(fac2, z, fac2);
    num_exp(fac2, fac2);
    num_mul(fac1, fac1, fac2);
    num_inv(_res, alpha);
    num_mul(fac1, fac1, _res);
    asymptotic_series(fac2, z, alpha, beta);
    num_add(_res, fac1, fac2);
    num_set(res, _res);
    delete(_res), delete(fac1), delete(fac2);
}

void
mittleff2 (num_t res,
           const num_t alpha, const num_t beta,
           const num_t z,
           const num_t acc)
{
    log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__,
              num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc));
    num_t _res;
    _res = new(num);
    num_set_d(_res, 0.0);
    asymptotic_series(_res, z, alpha, beta);
    num_set(res, _res);
    delete(_res);
}

/* Apply eq. (2.6) */
void
mittleff3_4 (num_t res,
             const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t acc,
             const int flag)
{
    num_t c, one, J, pi, aux, th, fac, fac1, fac2, fac3, fac4;

    c = new(num), one = new(num);
    pi = new(num), J = new(num), aux = new(num);
    fac = new(num), fac3 = new(num), fac4 = new(num);
    th = new(num), fac1 = new(num), fac2 = new(num);

    num_one(one);
    num_set_d(pi, M_PI);
    num_set_d_d(J, 0.0, 1.0);

    /* compute c(theta) */
    num_inv(th, alpha);
    num_pow(th, z, th);
    num_arg(th, th);
    num_sub(th, th, pi);
    num_mul(aux, J, th);
    num_exp(c, aux);
    num_sub(c, aux, c);
    num_add(c, c, one);
    num_mul_d(c, c, 2.0);
    num_sqrt(c, c);

    if (flag == 4)
        num_neg(c, c);

    num_sub(aux, one, beta);
    num_div(aux, aux, alpha);
    num_pow(fac1, z, aux);

    num_inv(aux, alpha);
    num_pow(aux, z, aux);
    num_exp(fac2, aux);

    num_inv(aux, alpha);
    num_abs(fac4, z);
    num_pow(aux, fac4, aux);
    num_mul_d(fac4, aux, 0.5);

    num_sqrt(aux, fac4);
    num_mul(aux, c, aux);
    num_erfc(fac3, aux);

    num_mul(fac, fac1, fac2);
    num_mul(fac, fac, fac3);
    num_div(fac, fac, alpha);
    num_mul_d(fac, fac, 0.5);

    asymptotic_series(aux, z, alpha, beta);

    num_add(res, fac, aux);

    delete(c), delete(one);
    delete(pi), delete(J), delete(aux);
    delete(fac), delete(fac3), delete(fac4);
    delete(th), delete(fac1), delete(fac2);
}

void
mittleff3 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    mittleff3_4(res, alpha, beta, z, acc, 3);
}

void
mittleff4 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    mittleff3_4(res, alpha, beta, z, acc, 4);
}

/* apply eqs. (4.25) and (4.26) */
void
mittleff5 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    num_t zero, rmax, aux, fac1, fac2, fac3, two, one, d, phi;

    one = new(num), two = new(num);
    rmax = new(num), aux = new(num);
    fac1 = new(num), fac2 = new(num), fac3 = new(num);
    d = new(num);
    zero = new(num);
    num_zero(zero);
    phi = new(num);
    
    num_set_d(two, 2.0);
    num_one(one);
    
    if (num_le_d(beta, 1.0))
    {
        /* Compute rmax */
        num_zero(rmax);
        if (num_ge_d(beta, 0.0))
        {
            num_abs(fac1, z);
            num_mul(fac1, fac1, two);
            
            num_pow(fac2, two, alpha);

            num_pow(aux, two, beta);
            num_mul(aux, aux, acc);
            num_mul_d(aux, aux, M_PI/12.0);
            num_log(aux, aux);
            num_mul_d(aux, aux, -2.0);
            num_pow(fac3, aux, alpha);

            num_max3(rmax, fac1, fac2, fac3);
        }
        else
        {
            num_abs(fac1, beta);
            num_add(fac1, fac1, one);
            num_mul(fac1, fac1, two);
            num_pow(fac1, fac1, alpha);

            num_abs(fac2, z);
            num_mul(fac2, fac2, two);

            num_pow(aux, two, beta);
            num_mul(aux, aux, acc);

            num_abs(aux, beta);
            num_mul_d(d, aux, 4.0);
            num_pow(d, d, aux);
            num_add(aux, aux, two);
            num_mul(d, d, aux);
            num_mul_d(d, d, 12.0);
            
            num_div(aux, aux, d);
            num_log(aux, aux);
            num_mul_d(aux, aux, -2.0);
            num_pow(fac3, aux, alpha);

            num_max3(rmax, fac1, fac2, fac3);
        }

        num_t from, to;
        from = new(num), to = new(num);
        num_set_d(from, 0.0);
        num_set(to, rmax);
        A(fac1, z, alpha, beta, zero);
        num_mul_d(phi, alpha, M_PI);
        integrate_B(fac2, alpha, beta, z, phi, from, to);
        num_add(res, fac1, fac2);
        delete(from), delete(to);
    }
    else
    {
        if (num_ge_d(beta, 0.0))
        {
        }
        else
        {
        }

        num_t from, to;
        from = new(num), to = new(num);
        num_set_d(from, 0.5);
        num_set_d(to, 2.0 * num_to_d(rmax));
        A(fac1, z, alpha, beta, zero);
        num_mul_d(phi, alpha, M_PI);
        integrate_B(fac2, alpha, beta, z, phi, from, to);
        num_add(res, fac1, fac2);
        delete(from), delete(to);
        
    }
    delete(zero);
    delete(aux);
    delete(one), delete(two), delete(d);
    delete(rmax);
    delete(fac1), delete(fac2), delete(fac3);
    delete(phi);
    /* const double x = num_real_d(_z); */
    /* const double y = num_imag_d(_z); */
    /* const double alpha = num_to_d(_alpha); */
    /* const double beta = num_to_d(_beta); */
    /* const double eps = num_to_d(_acc); */

    /* /\* Compute r_max, equation (4.53) *\/ */
    /* double rmax = 0.0; */
    /* double abs_z = sqrt(x*x + y*y); */
    /* if (beta >= 0.0) */
    /*     rmax = MAX3(1.0, */
    /*                 pow(2.0, abs_z), */
    /*                 pow(-log(M_PI * eps/6.0), alpha)); */
    /* else */
    /*     rmax = MAX3(pow(fabs(beta) + 1.0, alpha), */
    /*                 2.0*abs_z, */
    /*                 pow(-2.0 * log(M_PI*eps/(6.0*(fabs(beta)+2.0)*pow(2.0*fabs(beta), beta))), alpha)); */

    /* num_t _res = new(num), a = new(num); */
    /* A(a, x, y, alpha, beta, 0.0); */
    /* if (beta <= 1.0) /\* Equation (4.25) *\/ */
    /* { */
    /*     num_t integ_b = new(num); */
    /*     integrate_B(integ_b, alpha, beta, x, y, M_PI * alpha, 0.0, rmax); */
    /*     num_add(_res, a, integ_b); */
    /*     delete(integ_b); */
    /* } */
    /* else /\* Equation (4.26) *\/ */
    /* { */
    /*     num_t integ_b = new(num), integ_c = new(num); */
    /*     integrate_B(integ_b, alpha, beta, x, y, M_PI * alpha, 0.5, rmax); */
    /*     integrate_C(integ_c, alpha, beta, x, y, 0.5, -M_PI * alpha, M_PI * alpha); */
    /*     num_add(_res, integ_b, integ_c); */
    /*     num_add(_res, _res, a); */
    /*     delete(integ_b), delete(integ_c); */
    /* } */
    /* num_set(res, _res); */
    /* delete(_res), delete(a); */
}

void
mittleff6 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    /* const double eps = 0.5; */
    /* /\* Compute r_max, equation (4.54) *\/ */
    /* double rmax = 0.0; */
    /* double abs_z = sqrt(x*x + y*y); */
    /* if (ge(beta, 0.0)) */
    /* { */
    /*     rmax = MAX3(pow(2.0, alpha), */
    /*                 2.0 * abs_z, */
    /*                 pow(-log(M_PI*eps*pow(2.0, beta)/12.0), alpha)); */
    /* } */
    /* else // TODO Check with Hansjoerg whether the brackets mean "ceil" here */
    /*     rmax = MAX3( */
    /*         pow(ceil(2.0 * (fabs(beta) + 1.0)), alpha), // HERE */
    /*         2.0 * abs_z, */
    /*         pow(ceil(-log((M_PI * pow(2.0, beta)*eps)/(12*(fabs(beta) + 2)*pow(4 * fabs(beta), fabs(beta))))), alpha)); */

    /* num_t _res = new(num); */
    /* if (le(beta, 1.0)) /\* Equation (4.31) *\/ */
    /* { */
    /*     num_t integ_b = new(num); */
    /*     integrate_B(integ_b, alpha, beta, x, y, 2.0 * M_PI * alpha/3.0, 0.0, rmax); */
    /*     num_set(_res, integ_b); */
    /*     delete(integ_b); */
    /* } */
    /* else /\* Equation (4.32) *\/ */
    /* { */
    /*     num_t integ_b = new(num), integ_c = new(num); */
    /*     integrate_B(integ_b, alpha, beta, x, y, 2.0 * M_PI * alpha/3.0, 0.5, rmax); */
    /*     integrate_C(integ_c, alpha, beta, x, y, 0.5, -2.0 * M_PI * alpha/3.0, 2.0 * M_PI * alpha/3.0); */
    /*     num_add(_res, integ_b, integ_c); */
    /*     delete(integ_b), delete(integ_c); */
    /* } */
    /* num_to_d_d(res, _res); */
    /* delete(_res); */
}


/* rhs of equation (2.3) */
void
asymptotic_series (num_t res, const num_t z, const num_t alpha, const num_t beta)
{
    int k, kmax;
    double abs_z;
    num_t absz, sum, tmp, fac1, fac2, one;

    absz = new(num), one = new(num);
    sum = new(num), tmp = new(num);
    fac1 = new(num), fac2 = new(num);

    num_one(one);
    num_abs(absz, z);

    /* compute kmax */
    num_inv(tmp, alpha);
    num_pow(fac1, absz, tmp);
    num_mul(fac1, tmp, fac1);
    num_ceil(fac1, fac1);
    num_add(fac1, fac1, one);    
    kmax = (int) num_to_d(fac1); //(ceil((1.0/alpha) * pow(abs_z, 1.0/alpha)) + 1.0);
    //log_trace("[%s] |z|=%.5e, kmax=%d", __func__, num_to_d(absz), kmax);

    /* -sum([z**(-k) * mp.rgamma(beta - alpha*k) for k in range(1, kmax + 1)]) */
    for (k = 1; k <= kmax; k++)
    {
        //log_trace("[%s] k=%d", __func__, k);
        /* fac1 <- z**(-k) */
        num_pow_d(fac1, z, (double) (-k));

        /* fac2 <- rgamma(beta - alpha * k) */
        num_mul_d(fac2, alpha, (double) k);
        num_sub(fac2, beta, fac2);
        num_rgamma(fac2, fac2);

        /* partial sum = fac1 * fac2 */
        num_mul(tmp, fac1, fac2);
        
        num_add(sum, sum, tmp);
    }
    num_set(res, sum);
    num_mul_d(res, res, -1.0);
    
    delete(absz), delete(one);
    delete(sum), delete(tmp);
    delete(fac1), delete(fac2);
}
