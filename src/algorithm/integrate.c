#include "integrate.h"

#include "log.h"
#include "new.h"

#include "flint/profiler.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "acb_modular.h"
#include "acb_calc.h"

int
f_sin(acb_ptr res, const acb_t z, void * param, slong order, slong prec);

num_t
integrate_B (void)
{

    acb_t s, t, a, b;
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = 64;
    goal = 0;
    
    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);
    mag_init(tol);

    acb_set_d(a, 0);
    acb_set_d(b, 100);
    acb_calc_integrate(s, f_sin, NULL, a, b, goal, tol, options, prec);

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tol);

    flint_cleanup_master();
                
    return new(num, 0.0, 0.0);
}

/* f(z) = sin(z) */
int
f_sin(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sin(res, z, prec);

    return 0;
}
