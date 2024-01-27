void
test_exp_0 (void)
{
    const double complex z = 0.8;
    const double complex expected = +2.22554092849247e+00+0.00000000000000e+00*I;
    const double complex computed = mittleff(1.0, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_exp_1 (void)
{
    const double complex z = 2.0;
    const double complex expected = +7.38905609893065e+00+0.00000000000000e+00*I;
    const double complex computed = mittleff(1.0, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_exp_2 (void)
{
    const double complex z = 3.0+4.0*I;
    const double complex expected = -1.31287830814622e+01-1.52007844630680e+01*I;
    const double complex computed = mittleff(1.0, 1.0, z);
    TEST_VALUE(expected, computed);
}
