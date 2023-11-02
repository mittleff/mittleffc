#include "algorithm.h"
#include "log.h"
#include "new.h"

num_t
mittleff0 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    const double _alpha = num_to_double(alpha);
    const double _beta = num_to_double(beta);
    const double _acc = num_to_double(acc);
    
     log_trace("[%s] (alpha, beta, z, acc) = (%g, %g, %g%+gj, %g)",
              __func__,
              _alpha,
              _beta,
              num_real_d(z),
              num_imag_d(z),
              _acc);
     num_t zero = new(num, 0.0, 0.0);
     num_t one = new(num, 1.0, 0.0);
     num_t two = new(num, 2.0, 0.0);

     const double abs_z = num_to_double(num_abs(z));
     const int k1 = (int) (ceil((2 - _beta)/_alpha) + 1);
     const int k2 = (int) (ceil(log(_acc*(1-abs_z))/log(abs_z)) + 1);
     const int kmax = (k1 > k2) ? k1 : k2;

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
