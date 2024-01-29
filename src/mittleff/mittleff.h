/*
 * This file is part of mittleffc (https://github.com/mittleff/mittleffc).
 *
 * mittleffc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * mittleffc is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * mittleffc. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file mittleff.h
 * @brief Interface for the main functions of the library.
 */
#ifndef __MITTLEFF_H__
#define __MITTLEFF_H__

/**
 * Computes the Mittag-Leffler function
 *
 * @brief Evaluation of the Mittag-Leffler (ML) function with 2 parameters.
 *
 * @details Use the algorithm developed in [1].  The routine evaluates an
 * approximation \f$E_t\f$ of the ML function \f$E\f$ such that 
 *
 * \f{equation}{
 * |E - E_t|/(1 + |E|) \approx \textrm{tol}
 * \f}.
 * \p alpha must be a real and positive scalar and \p beta a real scalar. The
 * ML function is defined as
 *
 * \f{equation}{
 * E_{\alpha,\beta}(z) = \sum_{k=0}^{\infty}\frac{z^k}{\Gamma(\alpha\,k+\beta)}
 * \f}
 *
 * REFERENCES
 *
 * [1] H. Seybold and R. Hilfer. *Numerical Algorithm for Calculating the
 * Generalized Mittag-Leffler Function*. SIAM Journal on Numerical Analysis,
 * 47(1), 69–88 (2008)
 *
 * @param[in,out] result Holds the final result of the computation
 * @param[in] a   Real positive escalar parameter.
 * @param[in] b    Real scalar parameter.
 * @param[in] x    Real part of the complex number argument.
 * @param[in] y    Imaginary part of the complex number argument.
 * @param[in] prec     Accuracy for the computation.
 */
int
mittleff_cmplx (double* res,
          const double a,
          const double b,
          const double x,
          const double y,
          const unsigned short int prec);

#endif /* __MITTLEFF_H__ */
