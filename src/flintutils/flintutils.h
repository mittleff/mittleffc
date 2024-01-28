/*
 * This file is part of mittleffc (https://github.com/mittleff/mittleffc).
 *
 * mittleffc is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * mittleffc is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * mittleffc. If not, see <https://www.gnu.org/licenses/>.
 */

/** 
 * @file flintutils.h
 * @brief Misc interface for interacting with FLINT.
 */
#ifndef __FLINTUTILS_H__
#define __FLINTUTILS_H__

#include <flint/arb.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944
#endif

/* Default precision for FLINT */
#define PREC 64

/**
 * Converts an arb_t object to double.
 */
double
arbtod (const arb_t x);

void
arb_fmod (arb_t res, const arb_t self, const arb_t other);

#endif /* __FLINTUTILS_H__ */
