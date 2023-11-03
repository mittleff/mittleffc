#include "sf.h"

#include <assert.h>
#include <arb.h>
#include <arb_hypgeom.h>

#include "arb2num.h"
#include "new.h"
#include "log.h"

num_t
num_erf (const num_t _self)
{
    assert(num_is_real(_self));

    num_t res;
    arb_t _res, self;

    arb_init(self); arb_init(_res);
    arb_set_d(self, num_to_double(_self));
    arb_hypgeom_erf(_res, self, 53);
    res = new(num, arbtod(_res), 0.0);
    arb_clear(self); arb_clear(_res);

    return res;
}
