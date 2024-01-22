void
test_z_zero_0 (void)
{
    const double complex expected = +1.00000000000000e+00;
    const double complex computed = mittleff(1.0, 1.0, 0.0);
    TEST_VALUE(expected, computed);
}

void
test_z_zero_1 (void)
{
    const double complex expected = +1.08912442105834e+00;
    const double complex computed = mittleff(1.0, 1.2, 0.0);
    TEST_VALUE(expected, computed);
}

void
test_z_zero_2 (void)
{
    const double complex expected = -9.35778720912873e-02;
    const double complex computed = mittleff(1.0, -0.1, 0.0);
    TEST_VALUE(expected, computed);
}
