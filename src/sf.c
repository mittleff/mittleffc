#include "sf.h"

#include <assert.h>

#include "new.h"

num_t
num_erf (const num_t x)
{
    assert(num_is_real(x));
    return new(num, 0.0, 0.0);
}
