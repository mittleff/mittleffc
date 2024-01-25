#ifndef _NUM_H_
#define _NUM_H_

/**
 * Defines a pointer to the ADT
 */
typedef struct num_struct_t * num_t;

num_t num_alloc (void);
void  num_free  (num_t self);

#endif /* _NUM_H_ */
