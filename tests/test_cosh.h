
void
test_cosh_0 (void)
{
    const double complex z = -1.42469226000000e+00-1.93094380000000e-01*I;
    const double complex expected = +2.15772016942189e+00+3.75733317193863e-01*I;
    const double complex computed = mittleff(2.0, 1.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_cosh_1 (void)
{
    const double complex z = +6.64679980000000e-01+1.46133978000000e+00*I;
    const double complex expected = 1.34270384074874e-01+7.10437973043830e-01*I;
    const double complex computed = mittleff(2.0, 1.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_cosh_2 (void)
{
    const double complex z = +8.94857660000000e-01+1.89827224000000e+00*I;
    const double complex expected = -4.59266292056142e-01+9.64999829685609e-01*I;
    const double complex computed = mittleff(2.0, 1.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_cosh_3 (void)
{
    const double complex z = -1.32190510000000e-01+2.22985315000000e+00*I;
    const double complex expected = -6.17729663704252e-01-1.04810500062338e-01*I;
    const double complex computed = mittleff(2.0, 1.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_cosh_4 (void)
{
    const double complex z = -1.60041154000000e+00+3.95133590000000e-01*I;
    const double complex expected = +2.37976084266865e+00-9.14839293944298e-01*I;
    const double complex computed = mittleff(2.0, 1.0, z*z);
    TEST_VALUE(expected, computed);
}
