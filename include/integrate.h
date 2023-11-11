#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

#include "num.h"

void
A (num_t res,
   const double rez,
   const double imz,
   const double alpha,
   const double beta,
   const double x);

void
B (num_t res,
   const double r,
   const double alpha,
   const double beta,
   const double x,
   const double y,
   const double phi);

void
integrate_B (num_t res,
             const double alpha,
             const double beta,
             const double x,
             const double y,
             const double phi,
             const double from,
             const double to);

void
integrate_C (num_t res,
             const double alpha,
             const double beta,
             const double x,
             const double y,
             const double rho,
             const double from,
             const double to);

#endif  /* __INTEGRATE_H__ */
