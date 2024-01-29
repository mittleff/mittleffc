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
 * @brief Implementation of interace with FLINT.
 */

#include <flintutils.h>
#include <flint/arf.h>
#include <math.h>

double
arbtod (const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

void
acbtod (double* res, const acb_t x)
{
    arb_t r, i;
    
    arb_init(r);
    arb_init(i);
    
    acb_get_real(r, x);
    acb_get_imag(i, x);

    res[0] = arbtod(r);
    res[1] = arbtod(i);

    arb_clear(r);
    acb_clear(i);
}

void
arb_fmod (arb_t res, const arb_t self, const arb_t other)
{
    const double _self = arbtod(self);
    const double _other = arbtod(other);
    
    arb_set_d(res, fmod(_self, _other));
}

double
acb_real_d(const acb_t x)
{
    double res;
    arb_t r;

    arb_init(r);    
    acb_get_real(r, x);
    res = arbtod(r);
    arb_clear(r);

    return res;
}

double
acb_imag_d(const acb_t x)
{
    double res;
    arb_t i;

    arb_init(i);    
    acb_get_imag(i, x);
    res = arbtod(i);
    arb_clear(i);

    return res;
}
