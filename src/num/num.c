#include <num.h>
#include <acb.h>

#include <assert.h>

/* abstract data type */
struct num_struct_t {
    acb_t dat;
};
    
num_t
num_alloc (void)
{
    num_t self = malloc(sizeof(struct num_struct_t));
    assert (self != NULL);
    acb_init (self -> dat);
    return self;
}

void
num_free (num_t self)
{
    acb_clear (self -> dat);
    free(self);
}
