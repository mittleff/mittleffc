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
 * @file quad.h
 * @brief Interface for the quadrature routines.
 */
#ifndef __QUAD_H__
#define __QUAD_H__

#include "num.h"

#include "flint/acb.h"
#include "flint/acb_calc.h"

typedef struct
{
    void (* function) (num_t res, const num_t x, void * params);
    void * params;
} num_function_t;

//typedef void * num_function_t;

/* int */
/* f_integrand (acb_ptr res, const acb_t z, void * params, slong order, slong prec); */


void
quad (num_t res,
      num_function_t * F,
      const num_t from,
      const num_t to);

#endif /* __QUAD_H__ */
