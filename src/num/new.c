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
 * @file new.c
 * @brief Implementation of the initialization routines of the classes.
 */
#include <assert.h>
#include <stdlib.h>

#include "new.h"
#include "abc.h"

void*
new (const void* _class, ...)
{
    const struct ABC* class = _class;
	void* p = calloc(1, class->size);

	assert(p);
	* (const struct ABC **) p = class;

	if (class->ctor)
	{
        va_list ap;

		va_start(ap, _class);
		p = class->ctor(p, & ap);
		va_end(ap);
	}
	return p;
}

void
delete (void* self)
{
    const struct ABC** cp = self;

	if (self && (*cp) && (*cp)->dtor)
		self = (*cp)->dtor(self);
    
	free(self);
}

void*
clone (const void* self)
{
    const struct ABC* const* cp = self;

	assert(self && (*cp) && (*cp)->clone);
    
	return (*cp)->clone(self);
}

size_t
size_of (const void* self)
{
    const struct ABC* const* cp = self;

	assert(self && (*cp));
    
	return (*cp) -> size;
}
