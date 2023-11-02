#include "partition.h"

#include "log.h"

#include "new.h"
#include "num.h"

#define R0 0.95

static bool diskp (const num_t z, const num_t r);
static bool closure_diskp (const num_t z, const num_t r);

bool
in_region_G0 (const num_t z)
{
    /* log_trace("[%s] called with parameters: z=(%+.5e, %+.5e)", */
    /*           __func__, num_real_d(z), num_imag_d(z)); */
    return closure_diskp(z, new(num, R0, 0.0));
}

static bool
diskp (const num_t z, const num_t r)
{
    return num_lt(num_abs(z), r);
}

static bool
closure_diskp (const num_t z, const num_t r)
{
    return diskp(z, r) || num_eq(num_abs(z), r);
}
