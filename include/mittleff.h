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
#ifndef __ML_MITTLEFF_H__
#define __ML_MITTLEFF_H__

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
 * @param[in,out] res Holds the final result of the computation
 * @param[in] alpha   Real positive escalar parameter.
 * @param[in] beta    Real scalar parameter.
 * @param[in] re_z    Real part of the complex number argument.
 * @param[in] im_z    Imaginary part of the complex number argument.
 * @param[in] tol     Accuracy for the computation.
 */
int
mittleff_cmplx (double* res,
                const double alpha, const double beta,
                const double re_z, const double im_z,
                const double tol);

#endif /* __ML_MITTLEFF_H__ */
