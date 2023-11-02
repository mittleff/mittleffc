void
test_Memory_Management (void)
{
    num_t x;

    x = new(num);
    
    TEST_ASSERT_NOT_NULL(x);

    delete(x);
}

void
test_Initialization_and_Accessors (void)
{
    num_t x, y;

    x = new(num);
    TEST_ASSERT_NOT_NULL(x);

    num_set_d (x, +2.0);
    TEST_ASSERT_EQUAL_DOUBLE( +2.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_imag_d (x) );

    num_set_d_d (x, +3.0, -4.0);
    TEST_ASSERT_EQUAL_DOUBLE( +3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (x) );

    y = new(num);
    num_set_num (y, x);
    TEST_ASSERT_EQUAL_DOUBLE( +3.0, num_real_d (y) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (y) );

    delete(x), delete(y);
}

void
test_Precision_and_Comparisons (void)
{
    num_t x, y;

    x = new(num), y = new(num);
    TEST_ASSERT_NOT_NULL(x);
    TEST_ASSERT_NOT_NULL(y);
    TEST_ASSERT(num_is_zero(x));
    TEST_ASSERT(num_is_zero(y));

    /* x = 3.0, y = 3.0, z = x */
    num_set_d (x, 3.0), num_set_d (y, 3.0);
    TEST_ASSERT(num_eq (x, y));

    num_set_d_d (x, 2.0, -10.0);
    TEST_ASSERT_FALSE(num_eq (x, y));

    num_set_d (x, 0.0), num_set_d (y, 0.0);
    TEST_ASSERT(num_eq (x, y));

    num_set_d (x, 2.0);
    TEST_ASSERT( num_is_real (x) );

    num_set_d_d (x, 2.0, -10.0);
    TEST_ASSERT_FALSE( num_is_real (x) );

    num_set_d (x, 2.0), num_set_d (y, -3.0);
    TEST_ASSERT_TRUE( num_gt (x, y) );
    TEST_ASSERT_TRUE( num_lt (y, x) );
    TEST_ASSERT_TRUE( num_ge (x, y) );
    TEST_ASSERT_TRUE( num_le (y, x) );

    num_set_d (x, 2.0), num_set_d (y, 2.0);
    TEST_ASSERT_FALSE( num_gt (x, y) );
    TEST_ASSERT_FALSE( num_lt (y, x) );
    TEST_ASSERT_TRUE( num_ge (x, y) );
    TEST_ASSERT_TRUE( num_le (y, x) );

    /* Issue */
    /* num_set_d_d (x, 4.94623e-05, +4.62743e-21); */
    /* TEST_ASSERT_TRUE( num_is_real (x) ); */

    delete(x), delete(y);
}

void
test_Arithmetic (void)
{
    num_t x, y, z;

    x = new(num), y = new(num);
    TEST_ASSERT_NOT_NULL(x);
    TEST_ASSERT_NOT_NULL(y);

    num_set_d_d (x, 3.0, -4.0);
    TEST_ASSERT_EQUAL_DOUBLE( +3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (x) );
    num_neg (x, x);
    TEST_ASSERT_EQUAL_DOUBLE( -3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( +4.0, num_imag_d (x) );
    num_conj (x, x);
    TEST_ASSERT_EQUAL_DOUBLE( -3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (x) );

    
    z = new(num);
    num_set_d_d (x, +3.0, +4.0), num_set_d_d (y, +3.0, -4.0);

    /* add */
    num_add (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE( +6.0, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_imag_d (z) );

    /* sub */
    num_sub (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE( +8.0, num_imag_d (z) );

    /* mul */
    num_mul (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE( +25.0, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE( 0.0, num_imag_d (z) );

    /* div */
    num_div (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE( -0.28, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE( +0.96, num_imag_d (z) );

    delete(x), delete(y), delete(z);
}

void
test_Special_Functions (void)
{
    num_t x, y;

    x = new(num), y = new(num);
    TEST_ASSERT_NOT_NULL(x);
    TEST_ASSERT_NOT_NULL(y);

    num_set_d (x, +1.0), num_rgamma (y, x);
    TEST_ASSERT_EQUAL_DOUBLE( +1.0, num_real_d (y) );
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_imag_d (y) );

    num_set_d (x, +4.0), num_rgamma (y, x);
    TEST_ASSERT_EQUAL_DOUBLE( 1.0/6.0, num_real_d (y) );
    TEST_ASSERT_EQUAL_DOUBLE( 0.0, num_imag_d (y) );

    delete(x), delete(y);
}
