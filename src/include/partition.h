#ifndef _PARTITION_H_
#define _PARTITION_H_

#include <complex.h>

unsigned short int
region (double complex z, double alpha, double taylor_radius, slong prec);

// Region checking functions for the Mittag-Leffler algorithm
// These functions determine which computational branch to use based on
// the complex argument z and parameter alpha

// C wrapper functions for use with standard C types
// These functions provide a simple interface for other languages

/* double compute_r1(double a, long int prec); */

/* /\** */
/*  * C wrapper for region G0 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G0, 0 otherwise */
/*  *\/ */
/* unsigned short int in_region_G0(double complex z, double alpha, double r0, long int prec); */

/* /\** */
/*  * C wrapper for region G1 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G1, 0 otherwise */
/*  *\/ */
/* unsigned short int in_region_G1(double complex z, double alpha, double r0, long int prec); */

/* /\** */
/*  * C wrapper for region G2 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G2, 0 otherwise */
/*  *\/ */
/* unsigned short int in_region_G2(double complex z, double alpha, double r0, long int prec); */

/* /\** */
/*  * C wrapper for region G3 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G3, 0 otherwise */
/*  *\/ */
/* unsigned short int in_region_G3(double complex z, double alpha, double r0, long int prec); */

/* /\** */
/*  * C wrapper for region G4 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G4, 0 otherwise */
/*  *\/ */
/* unsigned short int in_region_G4(double complex z, double alpha, double r0, long int prec); */

/* /\** */
/*  * C wrapper for region G5 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G5, 0 otherwise */
/*  *\/ */
/* unsigned short int in_region_G5(double complex z, double alpha, double r0, long int prec); */

/* /\** */
/*  * C wrapper for region G6 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G6, 0 otherwise */
/*  *\/ */
/* unsigned short int in_region_G6(double complex z, double alpha, double r0, long int prec); */

#endif /* _PARTITION_H_ */
