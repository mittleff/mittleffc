#include "algorithm.h"
#include "log.h"
#include "new.h"

num_t
mittleff0 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
     log_trace("[%s] (alpha, beta, z, acc) = (%g, %g, %g%+gj, %g)",
              __func__,
              num_to_double(alpha),
              num_to_double(beta),
              num_real_d(z),
              num_imag_d(z),
              acc);
     num_t zero = new(num, 0.0, 0.0);
     num_t one = new(num, 1.0, 0.0);
     num_t two = new(num, 2.0, 0.0);

     int kmax = (int) num_to_double(
         num_max(
             num_add(
                 one,
                 num_ceil(num_div(num_sub(two, beta), alpha))),
             num_add(
                 one,
                 num_ceil(
                     num_div(
                         num_log(num_mul(acc, num_sub(one, num_abs(z)))),
                         num_log(num_abs(z)))))));

     num_t sum = zero;
     for (int k = 0; k <= kmax; k++)
     {
         num_t k_num = new(num, (double)k, 0.0);
         sum = num_add(
             sum,
             num_mul(num_pow(z, k_num),
                     num_rgamma(num_add(num_mul(alpha, k_num), beta))));
     }

     delete(zero); delete(one); delete(two);
     
    return sum;
}
